import statistics
file = open("pr76_30.txt","r")
date = list()
for i in range(0,30):
    number = file.readline().split("\n")[0]
    date.append(float(number))

print(str(round(min(date),5)) + " & " + str(round(max(date),5))+ " & " + str(round(statistics.mean(date),5)) +  " & " +  str(round(statistics.median(date),5)) +  " & " + str(round(statistics.pstdev(date),5)))