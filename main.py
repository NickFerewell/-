import astropy
from astropy.io import ascii
import math


#ICRF2
icrf2File = open("icrf2-all.txt")
icrf2Table = icrf2File.readlines()
icrf2File.close()

ArcSecondsInSecond = 1 #15

icrf2DataStart = 23
icrf2Data = []

for i in range(icrf2DataStart, len(icrf2Table)):
    #Всё в угловных секундах
    RA = (60*60*int(icrf2Table[i][36:38]) + 60*int(icrf2Table[i][39:41]) + float(icrf2Table[i][42:53])) * ArcSecondsInSecond
    Dec = (-1 if icrf2Table[i][55] == "-" else 1)*( 60*60*int(icrf2Table[i][56:58]) + 60*int(icrf2Table[i][59:61]) + float(icrf2Table[i][62:72]) )

    dRA = float(icrf2Table[i][74:84]) * ArcSecondsInSecond
    dDec = float(icrf2Table[i][85:94])
    icrf2Data.append({
        "name": icrf2Table[i][5:21],
        "RA": RA,
        "Dec": Dec,
        "dRA": dRA,
        "dDec": dDec
    })

# print(icrf2Data[1])



#ICRF3
icrf3File = open("icrf3sx.txt")
icrf3Table = icrf3File.readlines()
icrf3File.close()

icrf3DataStart = 22
icrf3Data = []

for i in range(icrf3DataStart, len(icrf3Table)):
    #Всё в угловных секундах
    RA = (60*60*int(icrf3Table[i][40:42]) + 60*int(icrf3Table[i][43:45]) + float(icrf3Table[i][46:57])) * ArcSecondsInSecond
    Dec = (-1 if icrf3Table[i][61] == "-" else 1)*( 60*60*int(icrf3Table[i][62:64]) + 60*int(icrf3Table[i][65:67]) + float(icrf3Table[i][68:78]) )

    dRA = float(icrf3Table[i][83:93]) * ArcSecondsInSecond
    dDec = float(icrf3Table[i][98:107])
    icrf3Data.append({
        "name": icrf3Table[i][5:21],
        "RA": RA,
        "Dec": Dec,
        "dRA": dRA,
        "dDec": dDec
    })

# print(icrf3Data[0])


#Ищем общие объекты
commonObjects = []
for object2 in icrf2Data:
    for object3 in icrf3Data:
        if object2["name"] == object3["name"]:
            commonObjects.append({
                "name": object2["name"],
                "2RA": object2["RA"],
                "2Dec": object2["Dec"],
                "2dRA": object2["dRA"],
                "2dDec": object2["dDec"],
                "3RA": object3["RA"],
                "3Dec": object3["Dec"],
                "3dRA": object3["dRA"],
                "3dDec": object3["dDec"]
            })


dRAsum = 0
dDecsum = 0

for object in commonObjects:
    dRAsum += math.sqrt((object["3RA"] - object["2RA"])**2 + object["2dRA"]**2 + object["3dRA"]**2)
    dDecsum += math.sqrt((object["3Dec"] - object["2Dec"])**2 + object["2dDec"]**2 + object["3dDec"]**2)

dA = dRAsum / len(commonObjects)
dD = dDecsum / len(commonObjects)

print("dA = " + str(dA), "dD = " + str(dD), "угловых секунд")



for object in commonObjects:
    d0a = math.sqrt((object["3RA"] - object["2RA"])**2 + object["2dRA"]**2 + object["3dRA"]**2)
    d0d = math.sqrt((object["3Dec"] - object["2Dec"])**2 + object["2dDec"]**2 + object["3dDec"]**2)

    object["d1a"] = abs(d0a - dA)
    object["d1d"] = abs(d0d - dD)



#Деление на зоны по RA
objectsByRAZone = [[] for i in range(0,360)]
dRAbyRAZoneSum = [0] * 360
dRAbyRAZone = [0] * 360

dDecbyRAZoneSum = [0]*360
dDecbyRAZone = [0]*360

for object in sorted(commonObjects, key=lambda object: object["3RA"]):
    zone = math.floor(object["3RA"] / (60*60))
    objectsByRAZone[zone].append(object)
    dRAbyRAZoneSum[zone] += object["d1a"]
    dDecbyRAZoneSum[zone] += object["d1d"]

for i, zone in enumerate(objectsByRAZone):
    dRAbyRAZone[i] = dRAbyRAZoneSum[i] / (len(zone) if len(zone) != 0 else 1)
    dDecbyRAZone[i] = dDecbyRAZoneSum[i] / (len(zone) if len(zone) != 0 else 1)

for i, zone in enumerate(objectsByRAZone):
    for object in zone:
        object["d2a"] = abs(object["d1a"] - dRAbyRAZone[i])
        object["d2d"] = abs(object["d1d"] - dDecbyRAZone[i])



#Деление по зонам по Dec
objectsByDecZone = [[] for i in range(0,180)]
dRAbyDecZoneSum = [0] * 180
dRAbyDecZone = [0]*180

dDecbyDecZoneSum = [0]*180
dDecbyDecZone = [0]*180

for object in sorted(commonObjects, key=lambda object: object["3Dec"]):
    zone = math.floor(object["3Dec"] / (60*60) + 90)
    # print(object["3Dec"] / (60*60), zone)
    objectsByDecZone[zone].append(object)
    dRAbyDecZoneSum[zone] += object["d2a"]
    dDecbyDecZoneSum[zone] += object["d2d"]

for i, zone in enumerate(objectsByDecZone):
    dRAbyDecZone[i] = dRAbyDecZoneSum[i] / (len(zone) if len(zone) != 0 else 1)
    dDecbyDecZone[i] = dDecbyDecZoneSum[i] / (len(zone) if len(zone) != 0 else 1)

for i, zone in enumerate(objectsByDecZone):
    for object in zone:
        object["d3a"] = abs(object["d2a"] - dRAbyDecZone[i])
        object["d3d"] = abs(object["d2d"] - dDecbyDecZone[i])


SRA = 0
SDec = 0

for object in commonObjects:
    SRA += object["d3a"]
    SDec += object["d3d"]

SRA = SRA / len(commonObjects)
SDec = SDec / len(commonObjects)

print("ξα = " + str(SRA), "ξδ = " + str(SDec))


'''
dA = 0.00024977738415668727 секунд, dD = 0.004777604044945355 угловых секунд
ξα = 0.00031637340493857437 ξδ = 0.006921975971838758
'''