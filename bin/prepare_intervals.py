import sys

input_file = sys.argv[1]
interval_list = 'intervals.list'
interval_map = 'interval_map.txt'
ref_fa = 'ref.fa'
output_list = open(f'{interval_list}','w')
output_map = open(f'{interval_map}','w')
output_fa = open(f'{ref_fa}','w')

chr_count = 1

with open(input_file, 'r') as fh:
    for line in fh:
        line = line.strip()
        if '>' in line:
            header = line[1:]
            output_fa.write(f'>{chr_count}\n')
            output_list.write(f'{chr_count}\n')
            output_map.write(f'{chr_count}\t{header}\n')
            chr_count += 1
        else:
            output_fa.write(f'{line}\n')
output_fa.close()
output_map.close()
output_list.close()
