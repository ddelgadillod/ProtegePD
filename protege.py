#!/usr/bin/python3

# -*- coding: utf-8 -*-

"""
Created on Fri Jun 19 13:16:07 2020

@author: ddelgadillo
"""

from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
import phl as pt
import os
import sys
import platform
#from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO
import numpy as np

import dash
from dash import dcc
from dash import html
import time
import plotly.figure_factory as ff
from dash.dependencies import Input, Output
import logging
import warnings
import argparse

import subprocess
import shlex


warnings.filterwarnings("ignore")
#np.seterr(divide='ignore', invalid='ignore')
old_settings = np.seterr(all='ignore') 
logging.getLogger('werkzeug').setLevel(logging.ERROR)
cli = sys.modules['flask.cli']
cli.show_server_banner = lambda *x: None



parser = argparse.ArgumentParser(
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                description=''''                   PROTÉGÉ
                                
                                                PROTEin coding GEne for phylogenetic tag and identification
                                                Useful design and visualization of primers.
                                            '''
                                 )


parser.add_argument('-s','--seq', 
                    metavar='SEQ', 
                    type=str,
                    required=True,
                    dest = 'seqPath',
                    help='Fasta file withe genes sequences')

parser.add_argument('-c', '--consensus', 
                    metavar='CONS', 
                    type=float,
                    default = 90,
                    #required=False,
                    dest = 'consensusPerc',
                    help='Consensus percentage')

parser.add_argument('-g', '--nogapconsensus', 
                    #metavar='GAP', 
                    #required=False,
                    dest = 'gapConsensus',
                    help='No consider consensus with gaps',
                    action='store_false')

parser.add_argument('-d', '--codon', 
                    metavar='COD', 
                    type=int,
                    default = 7,
                    #required=False,
                    dest = 'codon',
                    help='Codon primer length')

parser.add_argument('-v','--verbose', 
                    help='increase output verbosity',
                    action='store_true')



args = parser.parse_args()

genFile = args.seqPath
consensusPerc =  args.consensusPerc
gapConsensus = args.gapConsensus
codon = args.codon
nucWindow = codon * 3

#genFile = 'filtered3_gyrB_genes.fas'
#genFile = '/supplementary/evv/PhyloTAGs_package/PhyloTAGs_package/test20_05/filtered3_gyrB_genes.fas'
#genFile = 'gyrB_genes.fas'


#print(args.accumulate(args.integers))

pPath, sPath = os.path.split(os.path.abspath(__file__))
sysInfo = platform.system()


#consensusPerc = 80
#gapConsensus = True
#codon = 7
#nucWindow = codon * 3

#genFile = 'filtered3_gyrB_genes.fas'
#genFile = '/supplementary/evv/PhyloTAGs_package/PhyloTAGs_package/test20_05/filtered3_gyrB_genes.fas'
#genFile = 'gyrB_genes.fas'



genFile = pPath + '/' + genFile
translatedName = 'translated_seqs_pL.fas'



idList = []
seqList = []
nucLen = []
aminoList = []
aminoLen = []
for seq_record in SeqIO.parse(genFile, 'fasta'):
    idList.append(seq_record.id)
    nuc = seq_record.seq._data.decode("utf-8")
    seqList.append(nuc)
    nucLen.append(len(nuc))
    nuc = Seq(nuc)
    amino = nuc.translate()
    amino = amino._data.decode("utf-8")
    aminoList.append(amino)
    aminoLen.append(len(amino))
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))



sequences = pd.DataFrame({'id' : idList,
                          'nuc_seq' : seqList,
                          'nuc_lenght' : nucLen,
                          'amino_seq' : aminoList,
                          'amino_lenght' : aminoLen})


sequences.to_csv('sequences.csv',sep = ';', index = True)

f = open(translatedName, 'w')
for i in range(0,len(sequences)):
    f.write('>' + sequences.id[i] + '\n')
    f.write(sequences.amino_seq[i] + '\n')
f.close()

muscle = pPath + '/muscle/muscle_lin'
in_file = translatedName
out_file =  'aligned_muscle_pl_'+ translatedName



#muscle_cline = MuscleCommandline(muscle, input=in_file, out=out_file)
muscle_cline  = str(muscle) + ' -in ' + str(in_file) + ' -out ' + out_file


print(muscle_cline)
#muscle_cline()



#muscle_cline = MuscleCommandline(muscle,input = in_file)
alnProc = subprocess.run([muscle, "-in", in_file, "-out", out_file], capture_output=True)
#stdout, stderr = subprocess.run([muscle, "-in", in_file, "-out", out_file], capture_output=True)
#align = AlignIO.read(StringIO(stdout), "fasta")
align = AlignIO.read(open(out_file), "fasta")
#print(align)


print("Alignment length %i" % align.get_alignment_length())
for record in align:
   print(record.seq + " " + record.id) 

idAl = []
alSqs = []
alSqsLen = []

for i in range(0, len(align._records)):
    idAl.append(align._records[i].id)
    alSqs.append(align._records[i].seq._data.decode('utf-8'))
    alSqsLen.append(len(align._records[i].seq._data))


alSequences = pd.DataFrame({'id' : idAl,
                          'al_amino_seq' : alSqs,
                          'al_amino_lenght' : alSqsLen})
                          

alSequences.to_csv('alSequences.csv',sep = ';', index = True)

sequences = pd.merge(sequences, alSequences, on='id')
aminoEq = []
for i in range(0, len(sequences)):
    tseq = sequences.nuc_seq[i]
    tamino = sequences.amino_seq[i]
    r = range(0,sequences.nuc_lenght[i] + 1 ,3)
    print(len(r))
    
    aminoCodon = []
    for j in range(0,sequences.amino_lenght[i]):
        ts = tseq[r[j]:r[j+1]]
        aminoCodon.append([tamino[j],ts])    
    aminoEq.append(aminoCodon)

AlbyCod = []
for j in range (0,len(sequences)):
    print(j)
    tEq = aminoEq[j]
    tAl = sequences.al_amino_seq[j]
    tAlbyCod = ''
    aminoCont = 0
    for i in range(0,len(tAl)):
        if tAl[i] == '-':
            tAlbyCod = tAlbyCod + '---'
        else:
            if tAl[i] == tEq[aminoCont][0]:
                tAlbyCod = tAlbyCod + str(tEq[aminoCont][1])
                aminoCont = aminoCont + 1
    AlbyCod.append(tAlbyCod)
    
sequences['AlbyCodons'] = AlbyCod


codAlLen = int(np.mean(sequences.al_amino_lenght*3))

genNumber = len(sequences)



sConsensus = ''
for i in range(0,codAlLen):
    pcList = [0,0,0,0,0]
    for h in range(0, len(AlbyCod)):       
        if AlbyCod[h][i] == '-':
            pcList[0] += 1
        elif AlbyCod[h][i] == 'G':
            pcList[1] += 1
        elif AlbyCod[h][i] == 'A':
            pcList[2] += 1
        elif AlbyCod[h][i] == 'C':
            pcList[3] += 1
        elif AlbyCod[h][i] == 'T':
            pcList[4] += 1
    sList = ['-','G','A','T','C']
    pcList = np.asarray(pcList)
    pcList = (pcList/genNumber) * 100
    pcList = pcList.round(2)
    orderList = list(np.argsort(-pcList))
    print('% ALIGNMENT POSITION ' + str(i) + ':')
    for j in range(0,len(sList)):
        print(str(sList[j]) + ' -> ' + str(pcList[j]) + '%')# + ', rank ' + str(orderList[j]))
    if not gapConsensus and pcList[0] > 0:
        sConsensus = sConsensus + '-'
        print('Consensus in position ' + str(i) + ' is -, asigned by gapConsensus = False')
    else:
        print('RANKING:')
        for j in range(0,len(orderList)):
            print(str(j+1) + '. ' + str(sList[orderList[j]]) + ' -> ' + str(pcList[orderList[j]])  + '%')
        if pcList[orderList[0]] >= consensusPerc:
            sConsensus = sConsensus + str(sList[orderList[0]])
            print('Consensus in position ' + str(i) + ' is ' + str(sList[orderList[0]]))
        elif pcList[0] >= (100 - consensusPerc):
            sConsensus = sConsensus + str(sList[orderList[0]])
            print('Consensus in position ' + str(i) + ' is ' + str(sList[orderList[0]]) + ', GAP by default, no posible degeneracies consensus')
        elif pcList[orderList[0]] < consensusPerc:
            sList.remove(sList[0])
            pcList = np.delete(pcList, [0])
            orderList = list(np.argsort(-pcList))
            acc = pcList[orderList[0]]
            accNuc = sList[orderList[0]]
            for k in range(1,len(orderList)):
                acc += pcList[orderList[k]]
                accNuc = accNuc + sList[orderList[k]]
                if acc > consensusPerc:
                    accNuc = sorted(list(accNuc))
                    accNuc = ''.join(map(str,accNuc))
                    accDeg = pt.degEquivalent(accNuc)
                    sConsensus = sConsensus + accDeg
                    print('  Consensus with accumulate sum ' + str(acc) + '% and degeneracie ' + str(accDeg) + '(' +  str(accNuc) + ')')        
                    break
                else:
                    print('   No consensus with Accumulate, Accumulated Degeneracie is ' + str(accNuc))        
        else:
            print('No consensus in position ' + str(i))

phyloDF = pd.DataFrame({'position':[],'degeneracies':[],'forwardPrimer':[],'reversePrimer':[]})


for i in range(0,len(sConsensus) - nucWindow + 1):
    
    pos = str(i+1) + '-' + str(nucWindow+i)
    frwdPrimer = str(sConsensus[i:nucWindow+i])
    frwd = Seq(frwdPrimer)
    rvrsPrimer = frwd.reverse_complement()
    rvrsPrimer = rvrsPrimer._data.decode('utf-8')
    if '-' in frwdPrimer:
        degs = 0
    else:
        primerInfo = pt.primerDeg(frwdPrimer)
        degs = primerInfo.primerNP()
    tempPhylo = pd.DataFrame({'position':[pos],'degeneracies':[degs],'forwardPrimer':[frwdPrimer],'reversePrimer':[rvrsPrimer]})
    #phyloDF = phyloDF.append(tempPhylo, ignore_index=True)
    

    phyloDF = pd.concat([phyloDF, tempPhylo], ignore_index=True)

    print(str(i+1) + '-' + str(nucWindow+i) + '->' + str(sConsensus[i:nucWindow+i]) + ', degs = ' + str(degs))
    

phyloDF.to_csv('protege_consensus.csv',sep = ';', index = True)
print('Consensus by alignment position:')
print(sConsensus)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']


app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


#nPrimers = 10

i = 0
TApprox = 1




crrntTime = time.strftime("%Y-%m-%d_%H%M")
dirName = 'results_' + crrntTime
#os.mkdir(dirName)

#phyloDF = pd.read_csv('phyloTAGsPrimers.csv', sep = ' ')
phyloDF = phyloDF[phyloDF.degeneracies != 0]

#phyloDF = phyloDF[phyloDF.degeneracies <= 10000]


print("###############################################")
print("                    CONTROL 1                  ")
print("###############################################")

minLim = min(phyloDF.degeneracies)
maxLim = max(phyloDF.degeneracies)
maxInd = max(phyloDF.index)
logminLim = np.log10(minLim)
logmaxLim = np.log10(maxLim) + 0.02

degMin = minLim
degMax = maxLim
degOptions = list(phyloDF.degeneracies.unique())
degOptions.sort()


print("###############################################")
print("                    CONTROL 2                  ")
print("###############################################")

phyloFilteredDF = pt.filterDF(phyloDF, degMin, degMax)


scat = pt.posDegScatter(phyloFilteredDF)

zoomScat = pt.zoomDegScatter(phyloFilteredDF)


m = {}
for i in range(0,int(np.log10(maxLim*100))):
    m.update({10**i : str(10**i)})
#scat.write_html('scatter_allDeg.html', auto_open=True)

#scatF = pt.posDegScatter(phyloFilteredDF)
#scatF.write_html('scatter_filtered.html', auto_open=True)

phyloFilteredDF = phyloFilteredDF.reset_index(drop = True)

print("###############################################")
print("                    CONTROL 3                  ")
print("###############################################")

tempDist = ff.create_distplot([[0,0], [0,0]],
                              ['Forward Primer','Reverse Primer'],
                              curve_type = 'normal', 
                              show_rug = False, 
                              bin_size = 2
                             )
                             
                             
                             
print("###############################################")
print("                    CONTROL 4                  ")
print("###############################################")


tempDist.update_layout(paper_bgcolor = 'rgba(255,255,255,0.75)'
                       )


tempDist = pt.tempDist(phyloFilteredDF, i, TApprox)
#tempDist.write_html('temp_dist.html', auto_open=True)


colors = {
    'background': '#ffffff',
    'text': '#000000',
    'background2' : '#ffffff',
    'text2': '#000000'
}




#opt = ddPrimerOptions(phyloFilteredDF)


#fig = go.Figure(data=[go.Scatter(x=[1, 2, 3], y=[4, 1, 2])])

def transform_value(value):
    return 10 ** value


app.layout = html.Div(#style={'backgroundColor': colors['background']}, 
                      style={
                            #'background-image': 'url("assets/dna1.png")',
                            #'background-repeat': 'no-repeat',
                            #'background-position': 'right top',
                            #'background-size': '150px 100px'
                            },
                        children=[
                            html.H2(children='Phylotag Design',
                                    style={
                                    'textAlign': 'center',
                                    'color': colors['text']
                                    }),
                            html.Div( 
                                    children=[
                                        html.H5(children='Primer degeneracies by position',
                                        style={
                                        'textAlign': 'center',
                                        'color': colors['text2']
                                        }),

                                    dcc.Graph(style={'backgroundColor': 'transparent'},
                                        id='scat',
                                        figure = scat,
                                    ),
                                    #html.Button('Download data CSV ', 
                                    #            id='submit-val',
                                    #            n_clicks=0),

                                    html.Div(id='output-container-range-slider'),                                    
                                    
                                    # dcc.RangeSlider(
                                        # id = 'phyloDF-range',
                                        # #step = None,
                                        # marks = [0, 100, 2000, 100000, 133406300],
                                        # min = logminLim,
                                        # max = logmaxLim,
                                        # step = 0.05,
                                        # value=[logminLim, logmaxLim],
                                        # #vertical = True
                                        # ),

                                  


                            ]),
                            
                            # html.Div([
                                 # dcc.RangeSlider(0, 3,
                                    # id='non-linear-range-slider',
                                    # marks={i: '{}'.format(10 ** i) for i in range(4)},
                                    # #value=[0.1, 2],
                                    # value=[logminLim, logmaxLim],
                                    # dots=False,
                                    # step=0.01,
                                    # updatemode='drag'
                                # ),
                            # html.Div(id='output-container-range-slider-non-linear', style={'marginTop': 20})
                            # ]),
                           
                                
                            html.Div([
                                
                                html.Pre(id='click-data'),
                                
                                dcc.Graph(
                                        id='zoomScat',
                                        figure = zoomScat,
                                    ),

                            html.Div(id='container-zoom-range-slider'),                                    
                    
                            # dcc.Slider(
                                # id = 'phyloDFZoom-range',
                                # min = 50,
                                # max = maxInd,
                                # step = 50,
                                # marks={i: ' ±{} bp'.format(i) for i in range(100, maxInd, 150)},
                                # value = int(maxInd/2)
                            # ),  

                            html.Pre(id='click-data-zoom'),

                            html.Pre(id = 'selected-primers'),
                                    #children='Phylotag Design',
                                    #style={
                                    #'textAlign': 'center',
                                    #'color': colors['text']
                                    #}),
                                
                            ]),
                            html.Div(
                            [
                                html.Button("Download CSV", id="btn_csv"),
                                dcc.Download(id="download-dataframe-csv"),
                            ]),
                            
                            html.Hr(),
                            html.Div(
                                    children=[
                                        
                                        dcc.RadioItems(
                                                    id = 'tApprox-radio',
                                                    options=[
                                                    {'label': 'Tm Wallace "Rule of thumb"', 'value': 1},
                                                    {'label': 'Approx 2 Based on GC content', 'value': 2},
                                                    {'label': 'Approx 3 Based on GC content', 'value': 3},
                                                    {'label': 'Nearest neighbor', 'value': 4}
                                                ],
                                                value = 1,
                                                style={'width': '50%', 'display': 'inline-block'}),

                                            #dcc.Dropdown(
                                                #id='primer-dropdown',
                                                #style={'width': '50%', 'display': 'inline-block'})
                                        
                                        
                            ]),
                            html.Div(id='tdd-output'),

                            #html.Div(id='pdd-output'),
                            html.Div(style={'backgroundColor': 'transparent'}, 
                                    children=[
                                        html.H5(children='Melting Temperature Distribution',
                                        style={
                                        'textAlign': 'center',
                                        'color': colors['text2']
                                        }),

                                    dcc.Graph(style={'backgroundColor': 'transparent'},
                                        id='tdist-plot',
                                        figure = tempDist,
                                    ),
                            ]),
                            html.Div(id='filter-opt'),
                            
        
    #html.Div(id='output-container-range-slider')


])#cierra layout

print("###############################################")
print("                    CONTROL 5                  ")
print("###############################################")

    
@app.callback(
    dash.dependencies.Output('scat', 'figure'),
    [dash.dependencies.Input('phyloDF-range', 'value')])
def update_output(value):
    degMin = value[0]
    degMin = int(10**degMin)
    degMax = value[1]
    degMax = int(10**degMax)
    #print('min: ', str(degMin))
    #print('max: ', str(degMax))
    phyloFilteredDF = pt.filterDF(phyloDF, degMin, degMax)
    phyloFilteredDF = phyloFilteredDF.reset_index(drop = True)
    scat = pt.posDegScatter(phyloFilteredDF)
    return scat

print("###############################################")
print("                    CONTROL 6                  ")
print("###############################################")

@app.callback(
    dash.dependencies.Output('output-container-range-slider', 'children'),
    [dash.dependencies.Input('phyloDF-range', 'value')])
def update_output2(value):
    value[0] = int(10**value[0])
    value[1] = int(10**value[1])
    return 'You have selected "{}"'.format(value)


print("###############################################")
print("                    CONTROL 7                  ")
print("###############################################")

@app.callback(
    dash.dependencies.Output('tdd-output', 'children'),
    [dash.dependencies.Input('tApprox-radio', 'value')])
def update_output_6(value):
    return 'You have selected "{}"'.format(value)



print("###############################################")
print("                    CONTROL 8                  ")
print("###############################################")

@app.callback(
    Output('filter-opt', 'children'),
    [Input('phyloDF-range', 'value'),
     Input('scat', 'clickData'),
     Input('tApprox-radio', 'value')])

def update_output_8(lim, clickData, tApprox):
    try:
        x = clickData['points'][0]['x']
        #y = clickData['points'][0]['y']
        x = str(x)
        optF =  'limits are ' + str(int(10**lim[0])) + ' and ' + str(int(10**lim[1])) + ', primer ' + x + ' with T Approx ' + str(tApprox)
    except TypeError:
        optF = ''
    return optF


print("###############################################")
print("                    CONTROL 9                  ")
print("###############################################")

@app.callback(
    Output('click-data', 'children'),
    [Input('scat', 'clickData')])
def display_click_data(clickData):
    try:
        x = clickData['points'][0]['x']
        y = clickData['points'][0]['y']
        #print('X : ' + str(x) + ', Y :' + str(y))
        slctdPrimer = 'You choose primer on ' + str(x) + ' with ' + str(y) + ' degeneracies'
    except TypeError:
        #print('PRIMER NOT SELECTED')
        slctdPrimer = ''    
        #raise
    return slctdPrimer

@app.callback(
    Output('zoomScat', 'figure'),
    [Input('scat', 'clickData'),
     Input('phyloDFZoom-range', 'value')])
def display_click_data_2(clickData, zoomLimit):
    #zoomLimit = 300
    try:
        x = clickData['points'][0]['x']
        #print(type(zoomLimit))
        x = str(x)
        ind = phyloDF[phyloDF.position == x].index.tolist()
        ind = int(ind[0])
        
        if zoomLimit > ind:
            fInf = 0
        else:
            fInf = ind - zoomLimit
        
        if zoomLimit > (maxInd -ind):
            fSup = maxInd
        else:
            fSup = ind + zoomLimit
    
        zoomphyloDF = pt.filterDFbyIndex(phyloDF, fInf, fSup)        
        zoomScat = pt.zoomDegScatter(zoomphyloDF)
    except TypeError:
        zoomScat = pt.zoomDegScatter(phyloDF)
    return zoomScat


@app.callback(
    Output('click-data-zoom', 'children'),
    [Input('zoomScat', 'clickData')])
def display_click_data_4(clickData):
    try:
        x = clickData['points'][0]['x']
        y = clickData['points'][0]['y']
        txt = clickData['points'][0]['text']
        #print('X : ' + str(x) + ', Y :' + str(y))
        #return clickData
        #cData = json.loads(clickData)
        #print(type(cData))
        #print(cData)
        slctdPrimer = 'You choose reverse primer on ' + str(x) + ' with ' + str(y) + ' degeneracies, primer ' + str(txt)
    except:
        slctdPrimer = ''
    return slctdPrimer


@app.callback(
    dash.dependencies.Output('container-zoom-range-slider', 'children'),
    [dash.dependencies.Input('phyloDFZoom-range', 'value')])
def update_output_9(value):
    zoomMssg = 'Zoom on reverse primers ±'  + str(value) + ' bp.' 
    return zoomMssg

@app.callback(
    Output('selected-primers', 'children'),
    [Input('scat', 'clickData'),
     Input('zoomScat', 'clickData')])
def update_output_10(clickForward, clickReverse):
    try:
        frwrdPrimer = clickForward['points'][0]['text']
        #fPosition = clickForward['points'][0]['x']
        #fDegeneracies  = clickForward['points'][0]['y']
        rvrsPrimer = clickReverse['points'][0]['text']
        #rPosition = clickReverse['points'][0]['x']
        #rDegeneracies = clickReverse['points'][0]['y']
        #x = click['points'][0]['x']
        #y = clickData['points'][0]['y']
        #x = str(x)
        #ind = phyloDF[phyloDF.position == x].index.tolist()
        #ind = int(ind[0])
        #print(clickForward)
        #print(clickReverse)
        #Mssg = 'You choose '  + str(frwrdPrimer) + ' as forward primer in ' + str(fPosition) + ' with  ' + (fDegeneracies) + ' and ' + str(rvrsPrimer) + ' as reverse primer in ' + str(rPosition) + ' with  ' + (rDegeneracies)
        Mssg = 'You choose ' + str(frwrdPrimer) + ' --- ' + str(rvrsPrimer)
        
    except:
        Mssg = ''
    #print(Mssg)
    return Mssg


@app.callback(
    Output('tdist-plot', 'figure'),
     [Input('tApprox-radio', 'value'),
     Input('scat', 'clickData'),
     Input('zoomScat', 'clickData')])

def update_output_7(tApprox, clickForward, clickReverse):
    try:
        frwrdPrimer = clickForward['points'][0]['text']
        rvrsPrimer = clickReverse['points'][0]['text']
        frwrdPrimer = str(frwrdPrimer)
        rvrsPrimer = str(rvrsPrimer)
        
        frwrdPrimerInfo = pt.primerDeg(frwrdPrimer)
        rvrsPrimerInfo = pt.primerDeg(rvrsPrimer)
        
        if tApprox == 1:
            frwdTemp = frwrdPrimerInfo.TmWallace()
            rvrsTemp = rvrsPrimerInfo.TmWallace()
        elif tApprox == 2:
            frwdTemp = frwrdPrimerInfo.TmAp2()
            rvrsTemp = rvrsPrimerInfo.TmAp2()
        elif tApprox == 3:
            frwdTemp = frwrdPrimerInfo.TmAp3()
            rvrsTemp = rvrsPrimerInfo.TmAp3()
        elif tApprox == 4:
            frwdTemp = frwrdPrimerInfo.TmNN()
            rvrsTemp = rvrsPrimerInfo.TmNN()
        else:
            frwdTemp = frwrdPrimerInfo.TmWallace()
            rvrsTemp = rvrsPrimerInfo.TmWallace()
            
        tempDist = ff.create_distplot([frwdTemp, rvrsTemp],
                                      ['Forward Primer','Reverse Primer'],
                                      curve_type = 'normal', 
                                      show_rug = False, 
                                      bin_size = 2
                                     )
        tempDist.update_layout(paper_bgcolor = 'rgba(255,255,255,0.75)'
                               )
        
    except TypeError:
        tempDist = ff.create_distplot([[0,0], [0,0]],
                                      ['Forward Primer','Reverse Primer'],
                                      curve_type = 'normal', 
                                      show_rug = False, 
                                      bin_size = 2
                                      )
        tempDist.update_layout(paper_bgcolor = 'rgba(255,255,255,0.75)'
                               )

        #print('NO PRIMER CHOOSEN')    
    return tempDist
'''
if __name__ == '__main__':
    app.run_server(debug = False)
'''
@app.callback(
    Output('output-container-range-slider-non-linear', 'children'),
    Input('non-linear-range-slider', 'value'))
def update_output_11(value):
    transformed_value = [transform_value(v) for v in value]
    return 'Linear Value: {}, Log Value: [{:0.2f}, {:0.2f}]'.format(
        str(value),
        transformed_value[0],
        transformed_value[1]
    )


@app.callback(
    Output("download-dataframe-csv", "data"),
    Input("btn_csv", "n_clicks"),
    prevent_initial_call=True,
)
def func(n_clicks):
    return dcc.send_data_frame(phyloDF.to_csv, 'pd_protege_' + crrntTime + ".csv")

if __name__ == '__main__':
    app.run_server(host='0.0.0.0',debug=False,dev_tools_ui=False,dev_tools_props_check=False)
    
    
ind = phyloDF[phyloDF.position == '690-710'].index.tolist()



  

  