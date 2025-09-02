import csv
th = [500, 1400, 1500]
cnt_all = 0
cnt_th = {t: 0 for t in th}

with open("AV_orf.txt", newline='') as file:
    r = csv.DictReader(file, delimiter=',') 
    for row in r:
        coord = row['1B']
        if coord != 'NA-NA' and '-' in coord:
            try:
                start, end = map(int, coord.split('-'))
            except ValueError:
                continue
            length = end - start + 1
            cnt_all += 1
            for t in th:
                if length > t:
                    cnt_th[t] += 1

print(f"Количество генов 1B: {cnt_all}")
for t in sorted(th, reverse=True):
    print(f"Количество генов длиной > {t} nt: {cnt_th[t]}")
