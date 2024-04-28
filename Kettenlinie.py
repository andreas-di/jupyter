import math
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from stl import mesh


L = 136
x0 = 30
anzahlKloetze = 17
laengeDurchKlotzSimuliert = int(L/anzahlKloetze) #die Bogenlänge, die durch einen Klotz simuliert wird
breiteKloetze = 4


def a_bestimmen(unterGrenze,oberGrenze):
    a = (unterGrenze+oberGrenze)/2
    f_a = (L-2*a*math.sinh(1/a*x0))
    if( abs(f_a) < 1/10000000 ):
        return a
    elif (f_a < 0):
            return a_bestimmen(a,oberGrenze);
    else:
        return a_bestimmen(unterGrenze,a);
    
a = a_bestimmen(0,L)
#---
y0 = -a*math.cosh(1/a*x0)
#---
kettenlinie = lambda x: (a*np.cosh(1/a * x)+y0)
ableitung = lambda x: (np.sinh(1/a * x))
#---
x_werte_l = lambda l: (np.arcsinh(l/a) * a)
#---
def x_werte_holzstuecke():
    x_werte_holzstuecke = []
    for i in range (-int(L/2), int(L/2+1), laengeDurchKlotzSimuliert):
        x_werte_holzstuecke.append(x_werte_l(i))
    return x_werte_holzstuecke

x_werte = x_werte_holzstuecke()
#---
def y_werte_holzstuecke(xWerteListe):
    yWerte = []
    for i in xWerteListe:
        yWerte.append(kettenlinie(i))
    return yWerte

y_werte = y_werte_holzstuecke(x_werte)
#---
def sekantenSteigungsWinkel(p1, p2): #p1 und p2 werden als Parameter übergeben
    m = (kettenlinie(p1)-kettenlinie(p2))/(p1-p2) #Steigung wird bestimmt
    return (np.arctan(m)/(2 * math.pi) * 360) #bestimmung des Steigungswinkels und umrechnen von Radiant in Grad

def sekantenSteigungsWinkelListe(x_werte):
    steigungsWinkel = []
    steigungsWinkel.append( np.arctan( (kettenlinie(x_werte[0])-(-1)*kettenlinie(x_werte[1]))/(x_werte[0]-x_werte[1]) ) / (2*math.pi)*360 )
    for i in range(0,len(x_werte)-1): # von 0 bis len-1
        steigungsWinkel.append(sekantenSteigungsWinkel(x_werte[i], x_werte[i+1]))
    
    steigungsWinkel.append( np.arctan( (kettenlinie(x_werte[len(x_werte)-1])-(-1)*kettenlinie(x_werte[len(x_werte)-2]))/(x_werte[len(x_werte)-1]-x_werte[len(x_werte)-2]) ) / (2*math.pi)*360 ) ##
    return steigungsWinkel

sekantenSteigungswinkelWerteListe = sekantenSteigungsWinkelListe(x_werte)
#---
def winkelZwischenKloetzenHalbe(sekantenSteigungsWinkel):
    listeWinkel = []
    listeWinkel.append(sekantenSteigungsWinkel[0])
    for i in range(1,len(sekantenSteigungsWinkel)-2):
        listeWinkel.append((180 - sekantenSteigungsWinkel[i+1] + sekantenSteigungsWinkel[i])/2)
    listeWinkel.append(abs(sekantenSteigungsWinkel[len(sekantenSteigungsWinkel)-1]))
    return listeWinkel

winkelZwischenKloetzenHalbeWerte = winkelZwischenKloetzenHalbe(sekantenSteigungswinkelWerteListe)
#---
def sekantenLaengen(xWerteListe):
    sekantenLaengen = []
    for i in range(0, len(xWerteListe)-1):
        y1 = kettenlinie(xWerteListe[i])
        y2 = kettenlinie(xWerteListe[i+1])
        sekantenLaengen.append(math.sqrt((xWerteListe[i+1]-xWerteListe[i])**2 + (y2-y1)**2))
    return sekantenLaengen

sekantenLaengenWerte = sekantenLaengen(x_werte)
#---
def zusetzlicheLaengenstuecke(winkelInHolzStueck):
    langeSeite = []
    for i in range(0,len(winkelInHolzStueck)):
        langeSeite.append(2 / math.tan(math.radians(winkelInHolzStueck[i])))
    return langeSeite

zusetzlicheLaengenstueckeWerte = zusetzlicheLaengenstuecke(winkelZwischenKloetzenHalbeWerte)
#---
def laengeKloetzeAussenListe(sekantenLaengenListe, zusetzlicheLaengenListe):
    laengen = []
    for i in range(0,len(sekantenLaengenListe)):#1 bis len-1
        laengen.append(sekantenLaengenListe[i] + zusetzlicheLaengenListe[i] + zusetzlicheLaengenListe[i+1])
    return laengen

laengeKloetzeAussenWerteListe = laengeKloetzeAussenListe(sekantenLaengenWerte, zusetzlicheLaengenstueckeWerte)
#---
def laengeKloetzeInnenListe(sekantenLaengenListe, zusetzlicheLaengenListe):
    laengen = []
    for i in range(0,len(sekantenLaengenListe)):#1 bis len-1
        laengen.append(sekantenLaengenListe[i] - zusetzlicheLaengenListe[i] - zusetzlicheLaengenListe[i+1])
    return laengen

laengeKloetzeInnenWerteListe = laengeKloetzeInnenListe(sekantenLaengenWerte, zusetzlicheLaengenstueckeWerte)
#---
def funktionZeichnen(f, xmin, xmax, label=None):
    x_values = np.linspace(xmin, xmax, 300) #eine Liste an x-Werten wird erzeugt. 300 x-Werte im Bereich zwischen xmin und xmax
    y_values = f(x_values) #die zugehörigen y-Werte werden bestimmt

    plt.axhline(0, color='black',linewidth=1)
    plt.axvline(0, color='black',linewidth=1)

    plt.plot(x_values, y_values, label=label) #die Punkte, die die x- und y-Werte bilden werden eingezeichnet, da es relativ viele sind, sieht es nach einer Kurve aus

    plt.xlabel('x') #Achsenbeschriftung
    plt.ylabel('y')
    plt.title(f'Graph der Funktion {label}')
    plt.grid(True)
    plt.legend()
    plt.show()
#---
def kloetzeAbspeichern(sekantenLaengen, laengeKloetzeAussen, zusetzlicheLaenge, x_werte, y_werte, sekantenSteigungswinkelListe):
    for n in range (0, len(laengeKloetzeAussen)): #Schleife läuft von Null bis zur Anzahl der Klötze
        theta = np.radians(sekantenSteigungswinkelListe[n+1]) + np.radians(180) #theta ist der Winkel, um den der Klotz gedreht wird, bei n+1 weil trick mit erster und letzter klotz +180 da rechtsrum gedreht wird
        c, s = np.cos(theta), np.sin(theta) #in c wird der cosinus von theta geschrieben, in s der sin
        R = np.array(((c, -s), (s, c))) #Drehmatrix R wird initialisiert (siehe PDF Gleichung (7))

        #die Koordinaten der Punkte werden bestimmt (siehe PDF Abbildung 26). Hier liegt der Klotz noch an dem Nullpunkt
        a = np.array([-zusetzlicheLaenge[n+1], 2])
        b = np.array([ zusetzlicheLaenge[n+1],-2])
        c = np.array([-zusetzlicheLaenge[n+1]+laengeKloetzeAussen[n], 2])
        d = np.array([-zusetzlicheLaenge[n+1]+laengeKloetzeAussen[n]-2*zusetzlicheLaenge[n],-2])
        
        v = np.array([x_werte[n+1], y_werte[n+1]]) #v ist der Ortsvektor des Punktes, der beim bestimmen der Koordinaten auf den Nullpunkt gelegt wurde

        #die die Punkte a,b,c,d werden mit der Drehmatrix R multipliziert und mit dem Vektor v addiert und liegen jetzt auf der richtigen Stelle auf dem Funktionsgraphen
        a1 = np.dot(R,a)+v
        b1 = np.dot(R,b)+v
        c1 = np.dot(R,c)+v
        d1 = np.dot(R,d)+v

        plt.rcParams["figure.figsize"] = [7.00, 3.50] #legt die Größe der Matplotlib-Figuren fest (Größen in Zoll)
        plt.rcParams["figure.autolayout"] = True #aktivieren das automatische Layout

        polygon2 = Polygon([a1,c1,d1,b1,]) #eine Figur des Klotzes wird erstellt. Hierbei sind die Seiten a1c1 c1d1 d1b1 und b1a1

        x, y = polygon2.exterior.xy #der Variable x werden die x-Werte der Koordinaten zugeordnet, der Variable y die y-Werte
        plt.plot(x, y, c="black") #Die Punkte werden gezeichnet und gezeichnet (mit den Seiten des polygon2)

        #definition der Eckpunkte eines Klotzes a[0]/b[0]/... sind die x-Koordinaten der jeweiligen Punkte, a[1]/b[1]/... die y-Koordinaten
        #0 und vier sind die z-Koordinaten. Also ist ein Klotz hier 4LE hoch
        vertices = np.array([\
            [a[0], a[1] , 0], #[-1, -1, -1], #0
            [b[0], b[1] , 0], #[+1, -1, -1], #1
            [d[0], d[1] , 0], #[+1, +1, -1], #2
            [c[0], c[1] , 0], #[-1, +1, -1], #3
            [a[0], a[1] , 4], #[-1, -1, +1], #4 
            [b[0], b[1] , 4], #[+1, -1, +1], #5 
            [d[0], d[1] , 4], #[+1, +1, +1], #6
            [c[0], c[1] , 4]]) #[-1, +1, +1]]) #7
        #die Flächen des Klotzes werden durch Dreiecke definiert
        faces = np.array([\
            [0,3,1],
            [1,3,2],
            [0,4,7],
            [0,7,3],
            [4,5,6],
            [4,6,7],
            [5,1,2],
            [5,2,6],
            [2,3,6],
            [3,7,6],
            [0,1,5],
            [0,5,4]])

        cube = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype)) #ein Würfel mit den Flächen aus dem Array faces wird erzeugt
        for i, f in enumerate(faces): #schleife läuft über alle einträge im faces array
            for j in range(3): #schleife läuft über die drei Eckpunkte eines Dreiecks im faces Array
                cube.vectors[i][j] = vertices[f[j],:] #jetzt werden die Eckepunkte des cubes mit den Koordinaten aus dem vertices array gefüllt, also beschreibt dieser cube jetzt den aktuellen Klotz

        #der cube wird als stl datei gespeichert (dateityp der für den druck verwendet wird). Hier ist der Befehl auskommentiert, da sonst beim Ausführen der Funktion direkt ohen Nachfrage die stl Dateien der Klötze abgespeichert werden                
        #cube.save('cube'+str(n)+'.stl')

    #die Funktion der Kettenlinie wird gezeichnet (gleiches Verfahren wie in der Funktion zeichnen)
    x_values = np.linspace(-x0, x0, 300)
    y_values = kettenlinie(x_values)

    plt.axhline(0, color='black',linewidth=1)
    plt.axvline(0, color='black',linewidth=1)
    
    
    plt.plot(x_values, y_values)


    plt.show()
    
kloetzeAbspeichern(sekantenLaengenWerte, laengeKloetzeAussenWerteListe, zusetzlicheLaengenstueckeWerte, x_werte, y_werte, sekantenSteigungswinkelWerteListe)
