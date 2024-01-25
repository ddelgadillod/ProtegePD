#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 10:28:01 2020

@author: ddelgadillo
"""
import plotly.graph_objects as go
import plotly.figure_factory as ff
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq



class primerDeg:
    def __init__(self, primerDG):
        self.primerDG = primerDG
        
    def primerCheck(self):
        validPrimer = True
        for i in range(0,len(self.primerDG)):
            if not(self.primerDG[i] in ['G', 'A', 'T', 'C', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']):
                validPrimer = False
                break
        return validPrimer
    
    def primerNP(self):
        nc1 = ['G', 'A', 'T', 'C']
        nc2 = ['R', 'Y', 'M', 'K', 'S', 'W']
        nc3 = ['H', 'B', 'V', 'D']
        nc4 = ['N']
        pos = 0
        try:
            if self.primerCheck():
                pos = 1
                for i in range(0,len(self.primerDG)):
                    if self.primerDG[i] in nc1:
                        pos *= 1
                    elif self.primerDG[i] in nc2:
                        pos *= 2
                    elif self.primerDG[i] in nc3:
                        pos *= 3
                    elif self.primerDG[i] in nc4:
                        pos *= 4
            return pos
        except ValueError:
            print("Invalid Primer")
            return pos
    
    def primerND(self):
        nc2 = ['R', 'Y', 'M', 'K', 'S', 'W']
        nc3 = ['H', 'B', 'V', 'D']
        nc4 = ['N']
        try:
            if self.primerCheck():
                m = 0
                for i in range(0,len(self.primerDG)):
                    if (self.primerDG[i] in nc2) or (self.primerDG[i] in nc3) or (self.primerDG[i] in nc4):
                        m += 1
            return m
        except ValueError:
            print("Invalid Primer")
            return m
        
    def primerComb(self):
        nc1 = ['G', 'A', 'T', 'C']
        nc2 = ['R', 'Y', 'M', 'K', 'S', 'W']
        nc3 = ['H', 'B', 'V', 'D']
        nc4 = ['N']

        DGS = {'R':['G', 'A'], 'Y':['T', 'C'], 'M':['A', 'C'], 'K':['G', 'T'], 
               'S':['G', 'C'], 'W':['A', 'T'], 'H':['A', 'C', 'T'], 
               'B':['G', 'T', 'C'], 'V':['G', 'C', 'A'], 'D':['G', 'A', 'T'], 
               'N':['G', 'A', 'T', 'C']}

        LD = []
        LI = []
        LM = []
        lastIndex = self.primerNP()
        try:
            if self.primerCheck():
                for i in range(0,len(self.primerDG)):
                    if self.primerDG[i] in nc1:
                        LD.append(1)
                        LI.append(1)
                    else:
                        if self.primerDG[i] in nc2:
                            LD.append(2)
                            lastIndex= lastIndex//2
                        if self.primerDG[i] in nc3:
                            LD.append(3)
                            lastIndex= lastIndex//3
                        if self.primerDG[i] in nc4:
                            LD.append(4)
                            lastIndex= lastIndex//4
                        LI.append(lastIndex)
                    LM.append(self.primerNP()//(LD[i]*LI[i]))
            primerCombL = [None] * self.primerNP()
            for k in range(0,len(self.primerDG)):
                if LM[k] != self.primerNP():
                    G = []
                    opt = DGS[self.primerDG[k]]
                    for g in range(0,len(opt)):
                        for h in range(0,LI[k]):
                          G.append(opt[g])
                    G = G * LM[k]
                    for l in range(0,self.primerNP()):
                        if primerCombL[l] == None:
                            primerCombL[l] = G[k]
                        else:
                            primerCombL[l] = primerCombL[l] + G[l]
                else:
                    for l in range(0,self.primerNP()):
                        if primerCombL[l] == None:
                            primerCombL[l] = self.primerDG[k]
                        else:
                            primerCombL[l] = primerCombL[l] + self.primerDG[k]
            return primerCombL
        except ValueError:
            print("Invalid Primer")
            return primerCombL
    
    def primerBaseFreq(self):
        from collections import Counter as cntr
        baseFreq = []
        primerComb = self.primerComb()
        for i in range(0,len(primerComb)):
            degSeq = primerComb[i]
            degFreq = dict(cntr(degSeq))
            baseFreq.append(degFreq)
        return baseFreq
    
    def primerFreqG(self):
        primerBF = self.primerBaseFreq()
        pFreqGrouping = []
        for i in range(0,len(primerBF)):
            crrntPrimer = primerBF[i]
            G = 0
            C = 0
            A = 0
            T = 0
            if 'G' in crrntPrimer:
                G = crrntPrimer['G']
            if 'C' in crrntPrimer:
                C = crrntPrimer['C']
            if 'A' in crrntPrimer:
                A = crrntPrimer['A']
            if 'T' in crrntPrimer:
                T = crrntPrimer['T']
            pFreqGrouping.append({'GC':G + C, 'AT':A + T})
        return pFreqGrouping

    def GCporc(self):
        pFreqG = self.primerFreqG()
        GCporc = []
        for i in range(0,len(pFreqG)):
            crrntP = pFreqG[i]
            GCporc.append(round(crrntP['GC']/len(self.primerDG),5))
        return GCporc
    
    def TmWallace(self):
        TmW = []
        primerFreq = self.primerFreqG()
        for i in range(0,len(primerFreq)):
            crrntP = primerFreq[i]
            tTmW = 4 * crrntP['GC'] + 2 * crrntP['AT']
            TmW.append(tTmW)
        return TmW
            
    def TmAp2(self):
        Tm = []
        pFreqG = self.primerFreqG()
        for i in range (0,len(pFreqG)):
            crrntP = pFreqG[i]
            tTm = 64.9 + 41*((crrntP['GC'] - 16.4)/len(self.primerDG))
            #print('GC: ' + str(crrntP['GC']) + 'Tm: '  + str(tTm))
            Tm.append(round(tTm,4))
        return Tm
        
    def TmAp3(self):
        Tm = []
        primerComb = self.primerComb()
        for i in range(0,len(primerComb)):
            Tm.append(mt.Tm_GC(Seq(primerComb[i])))
        return Tm
        
    def TmNN(self):
        Tm = []
        primerComb = self.primerComb()
        for i in range(0,len(primerComb)):
            Tm.append(mt.Tm_NN(Seq(primerComb[i])))
        return Tm        
        
    
def posDegScatter(phyloDF):
    fig = go.Figure(data=go.Scatter(x = phyloDF['position'],
                                    y = phyloDF['degeneracies'],
                                    mode = 'markers',
                                    text = phyloDF['forwardPrimer'],
                                    #line_width = 1,
                                    line_color = 'rgb(34, 156, 0)',
                                    marker_size = 10,
                                    marker_opacity = 0.5
                                    )
                    )
    fig.update_layout(#width = 1000,
                      height = 530,
                      xaxis_title="Primer Position",
                      yaxis_title="Number of Degeneracies",
                      font_size = 11,
                      yaxis_type="log",
                      yaxis = dict( showexponent = 'all', 
                                   exponentformat = 'power'),
                      clickmode = 'event+select',
                      #opacity = 0.5,
                      #fig_bgcolor = "rgb(255, 255, 255)", 
                      #plot_bgcolor = "rgba(255, 255, 255, 255)", 
                      paper_bgcolor = 'rgba(255,255,255,0.75)'
                      #plot_bgcolor='rgba(0,0,0,0)',
                      #paper_bgcolor='rgba(0,0,0,0)'
                     )
   
    
    #fig.update_xaxes(tickangle=45)
    fig.update_xaxes(ticktext = phyloDF['position'])
    #fig.update_xaxes(showticklabels = False)

    #fig.show(config={'displayModeBar': True})
    #fig.show(config={'displayModeBar': False})
    return fig


def zoomDegScatter(phyloDF):
    fig = go.Figure(data=go.Scatter(x = phyloDF['position'],
                                    y = phyloDF['degeneracies'],
                                    text = phyloDF['reversePrimer'],
                                    mode = 'markers',
                                    #line_width = 1,
                                    line_color = 'rgb(255, 17, 0)',
                                    marker_size = 10,
                                    marker_opacity = 0.5
                                    )
                    )
    fig.update_layout(#width = 1000,
                      height = 530,
                      xaxis_title="Primer Position",
                      yaxis_title="Number of Degeneracies",
                      font_size = 11,
                      yaxis_type="log",
                      yaxis = dict( showexponent = 'all', 
                                   exponentformat = 'power'),
                      clickmode = 'event+select',
                      paper_bgcolor = 'rgba(255,255,255,0.75)'
                     )
   
    
    #fig.update_xaxes(tickangle=45)
    fig.update_xaxes(ticktext = phyloDF['position'])
    #fig.update_xaxes(showticklabels = False)

    #fig.show(config={'displayModeBar': True})
    #fig.show(config={'displayModeBar': False})
    return fig




def filterDF(DF,dMin,dMax):
    DF = DF[DF['degeneracies'] >= dMin]
    DF = DF[DF['degeneracies'] <= dMax]
    return DF

def filterDFbyIndex(DF,indMin,indMax):
    DF = DF[DF.index >= indMin]
    DF = DF[DF.index <= indMax]
    return DF

def primerInfo(DF, index, TApprox):
    primerStats = DF.loc[index]
    frwdInfoTemp = primerDeg(primerStats.forwardPrimer)
    rvrsInfoTemp = primerDeg(primerStats.reversePrimer)
    if TApprox == 1:
        frwdTemp = frwdInfoTemp.TmWallace()
        rvrsTemp = rvrsInfoTemp.TmWallace()
    else:
        frwdTemp = frwdInfoTemp.TmAp2()
        rvrsTemp = rvrsInfoTemp.TmAp2()

    TempDF = pd.DataFrame({'frwdTemp':frwdTemp,
                           'rvrsTemp':rvrsTemp})
    
    return TempDF


def ddPrimerOptions(DF):
    opt = []
    for ind in DF.index:
        lab = 'position ' + str(DF['position'][ind]) + ' Deg. ' + str(DF['degeneracies'][ind])
        val = ind
        tempDict = {'label': lab , 'value': val}
        opt.append(tempDict)
        #print(tempDict)
    return opt

def tempDist(DF, i, TApprox):
    TempDF = primerInfo(DF, i, TApprox)

    tempDist = ff.create_distplot([TempDF[c] for c in TempDF.columns], 
                                  ['Forward Dist', 'Reverse Dist'],
                                  #curve_type = 'normal', 
                                  show_rug = False, 
                                  bin_size = 2
                                 )
    tempDist.update_layout(
                      paper_bgcolor = 'rgba(255,255,255,0.75)'
                     )
    return tempDist


def degEquivalent(degeneracie):
    degDict = {'AG':'R',
               'CT':'Y',
               'CG':'S',
               'AT':'W',
               'GT':'K',
               'AC':'M',
               'CGT':'B',
               'AGT':'D',
               'ACT':'H',
               'ACG':'V',
               'ACGT':'N'}
    return degDict[degeneracie]
    