import csv
th = [500, 1400, 1500]
cnt_all = 0
cnt_th = {t: 0 for t in th}
with open("AV_orf.txt") as file:
    r = csv.DictReader(file, delimiter=',') 
    for row in r:
        coord = row['1B']
        if coord != 'NA-NA' and '-' in coord:
            try:
                start_end = coord.split('-')
                start = int(start_end[0])
                end = int(start_end[1])
            except ValueError:
                continue
            length = end - start + 1
            cnt_all += 1
            for t in th:
                if length > t:
                    cnt_th[t] += 1
print(f"Количество генов 1B: {cnt_all}")
for t in sorted(th):
    print(f"Количество генов длиной > {t} nt: {cnt_th[t]}")
