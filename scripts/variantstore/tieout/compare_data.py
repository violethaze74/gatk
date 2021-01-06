import sys

def get_next_line(i):
    for line in i:
        if (line.startswith("##")):
            pass;
        elif ("ZEROED_OUT_ASSAY" in line):
            pass;
        else:
            parts = line.strip().split("\t")
            loc = f"{parts[0]}:{parts[1]}"
            if (loc in exclude_list):
#                print(f"Skipping {loc}")
                pass;
            else:
                return line;

def compare_headers(header1, header2):
    a = set(header1.split("\t"))
    b = set(header2.split("\t"))
    diff = a.symmetric_difference(b)
    
    if (len(diff) == 0):
        print("Headers match, including all samples...")
    else:
        print(f"DIFF: headers are different! {a} vs {b}")
        sys.exit(1)

def parseline(e, header):
    data = {}
    
    parts = e.strip().split("\t")

    data['chrom'] = parts[0]
    data['pos'] = parts[1]
    data['id'] = parts[2]
    ref= parts[3]
    data['orig_alt'] = parts[4]
    
    alts = [x for x in parts[4].split(",")]
    
    ## and now minimize the ref and alt (with * allele being excempted)
    done = False
    while (not done and len(ref) != 1):
        if (all(ref[-1] == alt[-1] for alt in alts if alt != '*')):
            ref = ref[:-1]
            alts = [alt[:-1] if alt != '*' else '*' for alt in alts]
        else:
            done = True
    
    data['ref'] = ref
    data['alt'] = ",".join(alts)
    
    data['filter'] = parts[6]

    format_key = [x for x in parts[8].split(":")]

    samples = header.split("\t")[9:]
    sample_data = [ dict(zip(format_key, x.split(":"))) for x in parts[9:] ]
    
    data['sample_data'] =  dict(zip( samples, sample_data))

    return data;

def compare(e1, e2, key, ):
    if (e1[key] != e2[key]):
        print(f"DIFF on {key}")
        print(f"{e1}")
        print(f"{e2}")

def compare_float(e1, e2, key, tolerance):
    # compare directly first, also handles '.' case
    s1 = e1[key]
    s2 = e2[key]

    if (s1 != s2):
        if ("." in s1 or "." in s2):
            print(f"DIFF on {key} with values of {e1} and {e2}")

        else:
            v1 = float(s1)
            v2 = float(s2)

            delta = abs(v2 - v1)
            if delta > tolerance:
                print(f"DIFF on {key} of {delta}")
                print(f"{e1}")
                print(f"{e2}")

# need to sort, order doesn't matter (at this point)
def compare_alts(e1, e2):
    p1 = [x for x in e1.split(",") if x != '*']
    p1.sort()
    
    p2 = [x for x in e2.split(",") if x != '*']
    p2.sort()
    
    s1 = ",".join(p1)
    s2 = ",".join(p2)

    if (s1 != s2):
        print(f"DIFF on ALTS")
        print(f"{s1}")
        print(f"{s2}")

def get_gt_alleles(gt, ref, alts):
    alleles = [ref] + alts.split(",")
    delim = "|" if "|" in gt else "/"
    
    gt1 = gt.split(delim)[0]
    gt2 = gt.split(delim)[1]
        
    a1 = alleles[int(gt1)] if gt1 != "." else "."
    
#    if gt2 == "/":
#        print(gt)
        
    a2 = alleles[int(gt2)] if gt2 != "." else "."

    # TODO: for now, ignore phasing...
    return [a1,a2]
    
def compare_sample_data(e1, e2):
    sd1 = e1['sample_data']
    sd2 = e2['sample_data']
    
    if (len(sd1) != len(sd2)):
        print(f"DIFF on length of sample data {len(sd1)} and {len(sd2)}")
        print(f"{e1}")
        print(f"{e2}")
 
    for sample_id in sd1.keys():
        # if either has a GQ... compare it!
        if 'GQ' in sd1 or 'GQ' in sd2:
            compare(sd1[sample_id], sd2[sample_id], 'GQ')

        # compare genotypes based on actual alleles (since order of alts might differ)
        a1 = get_gt_alleles(sd1[sample_id]['GT'], e1['ref'], e1['alt'])
        a2 = get_gt_alleles(sd2[sample_id]['GT'], e2['ref'], e2['alt'])
        
        
        if (set(a1) != set(a2)):
            print(f"DIFF on Genotypes for {sample_id} at {e1['chrom']}:{e1['pos']} with {e1['alt']} and {e2['alt']}")
            print(f"{a1} vs {a2}")
            print(sd1[sample_id])
            print(sd2[sample_id])
            print("--------------")
    
    
# NOTE: files should have been passed through "unix sort" first
# i.e. cat foo.vcf | sort > new.vcf
vcf_file_1 = sys.argv[1]
vcf_file_2 = sys.argv[2]

exclude_list = []
if (len(sys.argv) == 4):
    with open(sys.argv[3]) as f:
        exclude_list = [x.strip() for x in f.readlines()]

lines = 0
with open(vcf_file_1) as file1, open(vcf_file_2) as file2:

    while True:
        line1 = get_next_line(file1)
        line2 = get_next_line(file2)

        if (line1 == None and line2 == None):
            break

        # get headers and compare (may have different sample order, that's ok)
        if (line1.startswith("#")):
            header1 = line1.strip()
            header2 = line2.strip()
            compare_headers(header1, header2)
            continue;
            
        # parse out data
        e1 = parseline(line1, header1)
        e2 = parseline(line2, header2)
        
        # do the comparison of exact matches
        for key in ['chrom','pos','id', 'ref']:
            compare(e1, e2, key)

        # TODO: temporary until we decide what to do with spanning deletions
        if ('*' in e1['alt']):
#            print(f"Dropping {e1['chrom']}:{e1['pos']} due to * allele")
            continue
            
        # compare the minimized version of ref/alt
        compare_alts(e1['alt'], e2['alt'])

        # compare sample level data
        compare_sample_data(e1, e2)
        
        lines = lines + 1
        
        
print(f"Compared {lines} lines...")
