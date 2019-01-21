# -*- coding: utf-8 -*-

from math import *
import numpy as np

class YXY():
    def __init__(self,gMat,gElt,gNdt,gBct):
        self.gMat = gMat
        self.gElt = gElt
        self.gNdt = gNdt
        self.gBct = gBct

    def DMtrElas(self,Mat):
        E=Mat[0]
        u=Mat[1]
        A1=E/(1-u*u)
        De=np.array([[A1,A1*u,0],[A1*u,A1,0],[0,0,A1*(1-u)/2]])
        return De

    def Tri3BMtr(self,xe):
        ai=[]
        bi=[]
        ci=[]
        xm=[]
        xm=xe
        xm.append(xe[0])
        xm.append(xe[1])
        for i in range(1,4):
            ai.append(xm[i][0]*xm[i+1][1]-xm[i+1][0]*xm[i][1])
            bi.append(xm[i][1]-xm[i+1][1])
            ci.append(xm[i+1][0]-xm[i][0])
        delta=sum(ai)
        A=delta/2
        if delta <= 1e-10:
            print('单元面积小于0')
        Bm=np.array([[bi[0],0,bi[1],0,bi[2],0],[0,ci[0],0,ci[1],0,ci[2]],[ci[0],bi[0],ci[1],bi[1],ci[2],bi[2]]])

        Bm=Bm/delta
        return Bm, A

    def Tri3EleStif(self,ie):
    #    MatNum=gEInfo[ie-1]
    #    Mat=gMat[MatNum-1] 如果gEInfo不是全是1怎么办
        De=self.DMtrElas(self.gMat)
        ENodes=[]
        ENodes=self.gElt[ie-1,:]
        xe=[]
        for i in range (0,len(ENodes)):
            xe.append(self.gNdt[int(ENodes[i])-1,:])
        [Bm,A]=self.Tri3BMtr(xe)
        Ke=np.dot(np.transpose(Bm),De)
        Ke=np.dot(Ke,Bm)
        Ke=np.dot(Ke,A)
        return Ke

    def AssembleEleStif(self,ie,Ke,gKA):
        ENodes=[]
        for i in range(0,len(self.gElt[ie-1])):
            ENodes.append(self.gElt[ie-1,i])
        enm=np.size(ENodes,0)
        dk=[]
        for i in range(0,enm):
            for j in range(0,enm):
                nodi=int(ENodes[i])
                nodj=int(ENodes[j])
                dk=gKA[2*nodi-2:2*nodi,2*nodj-2:2*nodj]
                gKA[2*nodi-2:2*nodi,2*nodj-2:2*nodj]=dk+Ke[2*i:2*i+2,2*j:2*j+2]       #坑！！！


    def LoopCalcGlobalStif(self,gKA):
        iLenElt=np.size(self.gElt,0)
        for ie in range(1,iLenElt+1):
            ENodes=np.zeros((len(self.gElt[ie-1,:]),1))
            ENodes[:,0]=self.gElt[ie-1,:]
            iN=np.size(ENodes,0)
            Ke=np.zeros((2*iN,2*iN))
            if iN==3:
                Ke=self.Tri3EleStif(ie)
            self.AssembleEleStif(ie,Ke,gKA)

    def LoopBoundaryDisp(self,gKA):
        iBct=np.size(self.gBct,0)
        for ib in range(1,iBct+1):
            iNode=int(self.gBct[ib-1,0])
            iType=int(self.gBct[ib-1,7])
            if iType==1:
                for ind in range(1,3):
                    irow=(iNode-1)*2+ind
                    if self.gBct[ib-1,ind]==1:
                        gKA[irow-1,:]=0.0
                        gKA[:,irow-1]=0.0
                        gKA[irow-1,irow-1]=1.0

    def LoopBoundaryLoad(self,gFA):
        iBct=np.size(self.gBct,0)
        for ib in range(1,iBct+1):
            iNode=int(self.gBct[ib-1,0])
            iType=int(self.gBct[ib-1,7])
            if iType==2:
                for ind in range(1,3):
                    irow=(iNode-1)*2+ind
                    gFA[irow-1]=gFA[irow-1]+self.gBct[ib-1,ind]
    def solve(self):
        gEInfo = []
        ielt = np.size(self.gElt, 0)
        for ie in range(0, ielt):
            gEInfo.append(1)
        iNdt = np.size(self.gNdt, 0)
        gKA = np.zeros((2 * iNdt, 2 * iNdt))
        gFA = np.zeros((2 * iNdt, 1))

        self.LoopCalcGlobalStif(gKA)
        self.LoopBoundaryLoad(gFA)
        self.LoopBoundaryDisp(gKA)

        Disp = np.dot(np.linalg.inv(gKA), gFA)
        gNTU = []
        for id in range(0, iNdt):
            gNTU.append(Disp[2 * id:2 * id + 2])
        gN = np.eye(22, 2)
        for i in range(0, 22):
            gN[i] = [gNTU[i][0][0], gNTU[i][1][0]]
        # print(gN)  # 变形值
        bianxing = self.gNdt + gN
        return bianxing
        # print(bianxing)  # 变形后的坐标
                   

import tkinter
from tkinter import *
import matplotlib.pyplot as plt


def getandplot():
    heig = int(text1.get())
    leng = int(text2.get())
    if leng/heig<8 or leng/heig>15:
        root1 = Tk()
        root1.geometry('100x100+0+0')
        Label(root1, text="请输入合适的梁长度").pack()
        root1.mainloop()
    else:
        huatu(heig,leng)

def huatu(h,l):
    le = l/10
    a = np.array([[0,0],
                  [0,h]])
    for i in range(20):
        if i%2==0:
            b = np.array([[(i/2+1)*le,0]])
            a = np.r_[a,b]
        else:
            b = np.array([[((i+1)/2)*le,h]])
            a = np.r_[a,b]

    plt.xlim(-0.2*l,1.2*l)
    plt.ylim(-0.4*l,0.3*l)
    for i in range(len(a) - 2):
        x = []
        x.append([a[i][0], a[i + 1][0]])
        x.append([a[i][0], a[i + 2][0]])
        y = []
        y.append([a[i][1], a[i + 1][1]])
        y.append([a[i][1], a[i + 2][1]])
        for j in range(2):
            plt.plot(x[j], y[j], color='r')
            plt.scatter(x[j], y[j], color='b')
    x = []
    x.append([a[-2][0], a[-1][0]])
    y = []
    y.append([a[-2][1], a[-1][1]])
    plt.plot(x[0], y[0], color='r')
    plt.show()

def jianzhi():
    f = int(text3.get())
    gBct = np.array([[1,1,1,0,0,0,0,1],
                     [21,1,1,0,0,0,0,1],
                     [12,0,-f,0,0,0,0,2]])
    gMat = np.array([800,0.3,0,0,0])
    h = int(text1.get())
    l = int(text2.get())
    le = l / 10
    gNdt = np.array([[0, 0],
                  [0, h]])
    for i in range(20):
        if i % 2 == 0:
            b = np.array([[(i / 2 + 1) * le, 0]])
            gNdt = np.r_[gNdt, b]
        else:
            b = np.array([[((i + 1) / 2) * le, h]])
            gNdt = np.r_[gNdt, b]
    gElt = np.array([[1, 3, 2],
                     [3, 4, 2],
                     [3, 5, 4],
                     [5, 6, 4],
                     [5, 7, 6],
                     [7, 8, 6],
                     [7, 9, 8],
                     [9, 10, 8],
                     [9, 11, 10],
                     [11, 12, 10],
                     [11, 13, 12],
                     [13, 14, 12],
                     [13, 15, 14],
                     [15, 16, 14],
                     [15, 17, 16],
                     [17, 18, 16],
                     [17, 19, 18],
                     [19, 20, 18],
                     [19, 21, 20],
                     [21, 22, 20]
                     ])

    yxy = YXY(gMat,gElt,gNdt,gBct)
    bianxing = yxy.solve()

    plt.xlim(-0.2*l,1.2*l)
    plt.ylim(-0.4*l, 0.3*l)
    for i in range(len(bianxing) - 2):
        x = []
        x.append([bianxing[i][0], bianxing[i + 1][0]])
        x.append([bianxing[i][0], bianxing[i + 2][0]])
        y = []
        y.append([bianxing[i][1], bianxing[i + 1][1]])
        y.append([bianxing[i][1], bianxing[i + 2][1]])
        for j in range(2):
            plt.plot(x[j], y[j], color='r')
            plt.scatter(x[j], y[j], color='b')
    x = []
    x.append([bianxing[-2][0], bianxing[-1][0]])
    y = []
    y.append([bianxing[-2][1], bianxing[-1][1]])
    plt.plot(x[0], y[0], color='r')

    arrowx = np.array([l/2, l/2])
    arrowy = np.array([1.2*h, 2 * h])
    arrowx1 = np.array([l/2, l/2 - 5])
    arrowy1 = np.array([1.2*h, 1.2*h + 5])
    arrowx2 = np.array([l/2, l/2 + 5])
    arrowy2 = np.array([1.2*h, 1.2*h + 5])
    plt.plot(arrowx, arrowy, color='r')
    plt.plot(arrowx1, arrowy1, color='r')
    plt.plot(arrowx2, arrowy2, color='r')
    plt.xlabel(" F = " + str(f) + "KN")
    plt.title("Simply supported beam")
    plt.show()

def xuanbi():
    f = int(text3.get())
    gBct = np.array([[1, 1, 1, 0, 0, 0, 0, 1],
                     [2, 1, 1, 0, 0, 0, 0, 1],
                     [22, 0, -f, 0, 0, 0, 0, 2]])
    gMat = np.array([1000, 0.25, 0, 0, 0])
    h = int(text1.get())
    l = int(text2.get())
    le = l / 10
    gNdt = np.array([[0, 0],
                     [0, h]])
    for i in range(20):
        if i % 2 == 0:
            b = np.array([[(i / 2 + 1) * le, 0]])
            gNdt = np.r_[gNdt, b]
        else:
            b = np.array([[((i + 1) / 2) * le, h]])
            gNdt = np.r_[gNdt, b]
    gElt = np.array([[1, 3, 2],
                     [3,4,2],
                     [3,5,4],
                     [5,6,4],
                     [5,7,6],
                     [7,8,6],
                     [7,9,8],
                     [9,10,8],
                     [9,11,10],
                     [11,12,10],
                     [11,13,12],
                     [13,14,12],
                     [13,15,14],
                     [15,16,14],
                     [15,17,16],
                     [17,18,16],
                     [17,19,18],
                     [19,20,18],
                     [19,21,20],
                     [21,22,20]
                     ])

    yxy = YXY(gMat, gElt, gNdt, gBct)
    bianxing = yxy.solve()

    plt.xlim(0,1.2*l)
    plt.ylim(-0.4*l, 0.3*l)
    for i in range(len(bianxing) - 2):
        x = []
        x.append([bianxing[i][0], bianxing[i + 1][0]])
        x.append([bianxing[i][0], bianxing[i + 2][0]])
        y = []
        y.append([bianxing[i][1], bianxing[i + 1][1]])
        y.append([bianxing[i][1], bianxing[i + 2][1]])
        for j in range(2):
            plt.plot(x[j], y[j], color='r')
            plt.scatter(x[j], y[j], color='b')
    x = []
    x.append([bianxing[-2][0], bianxing[-1][0]])
    y = []
    y.append([bianxing[-2][1], bianxing[-1][1]])
    plt.plot(x[0], y[0], color='r')

    arrowx = np.array([l,l])
    arrowy = np.array([h,1.8*h])
    arrowx1 = np.array([l,l-5])
    arrowy1 = np.array([h, h+5])
    arrowx2 = np.array([l,l+5])
    arrowy2 = np.array([h, h+5])
    plt.plot(arrowx,arrowy,color = 'r')
    plt.plot(arrowx1, arrowy1, color='r')
    plt.plot(arrowx2, arrowy2, color='r')
    plt.xlabel(" F = " + str(f) + "KN")
    plt.title("Cantilever")
    plt.show()


root = Tk()
root.geometry('600x400+0+0')
Label(root,text = '输入梁高(cm)，工程上常用20-50').pack()
text1 = Entry(root)
text1.pack()
Label(root,text = '输入梁长度(cm),建议为梁高的8-15倍').pack()
text2 = Entry(root)
text2.pack()
button1 =Button(root,text = '确定并显示',command = getandplot).pack()

Label(root,text = '输入荷载大小(KN)，建议值为50-120').pack()
text3 = Entry(root)
text3.pack()
button2  = Button(root,text = '生成简支梁（跨中荷载）变形图',command = jianzhi).pack()
button3 = Button(root,text = '生成悬臂梁（悬臂端荷载）变形图',command = xuanbi).pack()
root.mainloop()

