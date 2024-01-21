# -*- coding: utf-8 -*-

from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
import phl as pt
import os
import sys
import platform
from Bio.Align.Applications import MuscleCommandline
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

warnings.filterwarnings("ignore")
#np.seterr(divide='ignore', invalid='ignore')
old_settings = np.seterr(all='ignore')
logging.getLogger('werkzeug').setLevel(logging.ERROR)
cli = sys.modules['flask.cli']
cli.show_server_banner = lambda *x: None

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']


app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


#nPrimers = 10

i = 0
TApprox = 1

phyloDF = pd.read_csv('protege_consensus.csv', sep = ';')

phyloDF = phyloDF[phyloDF.degeneracies != 0]

minLim = min(phyloDF.degeneracies)
maxLim = max(phyloDF.degeneracies)
maxInd = max(phyloDF.index)
logminLim = np.log10(minLim)
logmaxLim = np.log10(maxLim) + 0.02

degMin = minLim
degMax = maxLim
degOptions = list(phyloDF.degeneracies.unique())
degOptions.sort()

print("CONTROL 2")

phyloFilteredDF = pt.filterDF(phyloDF, degMin, degMax)

scat = pt.posDegScatter(phyloFilteredDF)

zoomScat = pt.zoomDegScatter(phyloFilteredDF)

m = {}
for i in range(0, int(np.log10(maxLim * 100))):
    m.update({10 ** i: str(10 ** i)})
# scat.write_html('scatter_allDeg.html', auto_open=True)

# scatF = pt.posDegScatter(phyloFilteredDF)
# scatF.write_html('scatter_filtered.html', auto_open=True)

phyloFilteredDF = phyloFilteredDF.reset_index(drop=True)

##########################

tempDist = ff.create_distplot([[0, 0], [0, 0]],
                              ['Forward Primer', 'Reverse Primer'],
                              curve_type='normal',
                              show_rug=False,
                              bin_size=2
                              )
tempDist.update_layout(paper_bgcolor='rgba(255,255,255,0.75)'
                       )

print("CONTROL 3")

tempDist = pt.tempDist(phyloFilteredDF, i, TApprox)

tempDist.write_html('temp_dist.html', auto_open=True)


colors = {
    'background': '#ffffff',
    'text': '#000000',
    'background2': '#ffffff',
    'text2': '#000000'
}

# opt = ddPrimerOptions(phyloFilteredDF)


# fig = go.Figure(data=[go.Scatter(x=[1, 2, 3], y=[4, 1, 2])])


print("CONTROL 4")
app.layout = html.Div(  # style={'backgroundColor': colors['background']},
    style={
        # 'background-image': 'url("assets/dna1.png")',
        # 'background-repeat': 'no-repeat',
        # 'background-position': 'right top',
        # 'background-size': '150px 100px'
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
                          figure=scat,
                          ),
                # html.Button('Download data CSV ',
                #            id='submit-val',
                #            n_clicks=0),

                html.Div(id='output-container-range-slider'),

                dcc.RangeSlider(
                    id='phyloDF-range',
                    # step = None,
                    # marks = m,
                    min=logminLim,
                    max=logmaxLim,
                    step=0.05,
                    value=[logminLim, logmaxLim],
                    # vertical = True
                ),

            ]),

        html.Div([

            html.Pre(id='click-data'),

            dcc.Graph(
                id='zoomScat',
                figure=zoomScat,
            ),

            html.Div(id='container-zoom-range-slider'),

            dcc.Slider(
                id='phyloDFZoom-range',
                min=50,
                max=maxInd,
                step=50,
                marks={i: ' ±{} bp'.format(i) for i in range(100, maxInd, 150)},
                value=int(maxInd / 2)
            ),

            html.Pre(id='click-data-zoom'),

            html.Pre(id='selected-primers'),
            # children='Phylotag Design',
            # style={
            # 'textAlign': 'center',
            # 'color': colors['text']
            # }),

        ]),
        html.Hr(),
        html.Div(
            children=[

                dcc.RadioItems(
                    id='tApprox-radio',
                    options=[
                        {'label': 'Approx 1', 'value': 1},
                        {'label': 'Approx 2', 'value': 2}
                    ],
                    value=1,
                    style={'width': '50%', 'display': 'inline-block'}),

                # dcc.Dropdown(
                # id='primer-dropdown',
                # style={'width': '50%', 'display': 'inline-block'})

            ]),
        html.Div(id='tdd-output'),

        # html.Div(id='pdd-output'),
        html.Div(style={'backgroundColor': 'transparent'},
                 children=[
                     html.H5(children='Melting Temperature Distribution',
                             style={
                                 'textAlign': 'center',
                                 'color': colors['text2']
                             }),

                     dcc.Graph(style={'backgroundColor': 'transparent'},
                               id='tdist-plot',
                               figure=tempDist,
                               ),
                 ]),
        html.Div(id='filter-opt'),

        # html.Div(id='output-container-range-slider')

    ])  # cierra layout

print("CONTROL 5")


@app.callback(
    dash.dependencies.Output('scat', 'figure'),
    [dash.dependencies.Input('phyloDF-range', 'value')])
def update_output(value):
    degMin = value[0]
    degMin = int(10 ** degMin)
    degMax = value[1]
    degMax = int(10 ** degMax)
    # print('min: ', str(degMin))
    # print('max: ', str(degMax))
    phyloFilteredDF = pt.filterDF(phyloDF, degMin, degMax)
    phyloFilteredDF = phyloFilteredDF.reset_index(drop=True)
    scat = pt.posDegScatter(phyloFilteredDF)
    return scat


@app.callback(
    dash.dependencies.Output('output-container-range-slider', 'children'),
    [dash.dependencies.Input('phyloDF-range', 'value')])
def update_output2(value):
    value[0] = int(10 ** value[0])
    value[1] = int(10 ** value[1])
    return 'You have selected "{}"'.format(value)


@app.callback(
    dash.dependencies.Output('tdd-output', 'children'),
    [dash.dependencies.Input('tApprox-radio', 'value')])
def update_output_6(value):
    return 'You have selected "{}"'.format(value)


'''
@app.callback(
    Output('tdist-plot', 'figure'),
    [Input('phyloDF-range', 'value'),
     Input('scat', 'clickData'),
     Input('tApprox-radio', 'value')])

def update_output_7(lim, clickData, tApprox):
    degMin = lim[0]
    degMin = int(10**degMin)
    degMax = lim[1]
    degMax = int(10**degMax)
    #print('MIN: ', str(value[0]) + '-- MAX: ', str(value[1]))
    phyloFilteredDF3 = pt.filterDF(phyloDF, degMin, degMax)
    phyloFilteredDF3 = phyloFilteredDF3.reset_index(drop = True)
    x = clickData['points'][0]['x']
    #y = clickData['points'][0]['y']
    x = str(x)
    print('POSITION IS :' + x)

    primer = phyloDF[phyloFilteredDF3.position == x].index.tolist()
    primer = int(primer[0])
    print('PRIMER IS ' + str(primer))
    tempDist = pt.tempDist(phyloFilteredDF3, primer, tApprox)
    #optF = 'primer opt ' + str(primer) + ' with T Approx ' + str(tApprox)
    return tempDist
'''


@app.callback(
    Output('filter-opt', 'children'),
    [Input('phyloDF-range', 'value'),
     Input('scat', 'clickData'),
     Input('tApprox-radio', 'value')])
def update_output_8(lim, clickData, tApprox):
    try:
        x = clickData['points'][0]['x']
        # y = clickData['points'][0]['y']
        x = str(x)
        optF = 'limits are ' + str(int(10 ** lim[0])) + ' and ' + str(
            int(10 ** lim[1])) + ', primer ' + x + ' with T Approx ' + str(tApprox)
    except TypeError:
        optF = ''
    return optF


@app.callback(
    Output('click-data', 'children'),
    [Input('scat', 'clickData')])
def display_click_data(clickData):
    try:
        x = clickData['points'][0]['x']
        y = clickData['points'][0]['y']
        # print('X : ' + str(x) + ', Y :' + str(y))
        slctdPrimer = 'You choose primer on ' + str(x) + ' with ' + str(y) + ' degeneracies'
    except TypeError:
        # print('PRIMER NOT SELECTED')
        slctdPrimer = ''
        # raise
    return slctdPrimer


@app.callback(
    Output('zoomScat', 'figure'),
    [Input('scat', 'clickData'),
     Input('phyloDFZoom-range', 'value')])
def display_click_data_2(clickData, zoomLimit):
    # zoomLimit = 300
    try:
        x = clickData['points'][0]['x']
        # print(type(zoomLimit))
        x = str(x)
        ind = phyloDF[phyloDF.position == x].index.tolist()
        ind = int(ind[0])

        if zoomLimit > ind:
            fInf = 0
        else:
            fInf = ind - zoomLimit

        if zoomLimit > (maxInd - ind):
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
        # print('X : ' + str(x) + ', Y :' + str(y))
        # return clickData
        # cData = json.loads(clickData)
        # print(type(cData))
        # print(cData)
        slctdPrimer = 'You choose reverse primer on ' + str(x) + ' with ' + str(y) + ' degeneracies, primer ' + str(txt)
    except:
        slctdPrimer = ''
    return slctdPrimer


@app.callback(
    dash.dependencies.Output('container-zoom-range-slider', 'children'),
    [dash.dependencies.Input('phyloDFZoom-range', 'value')])
def update_output_9(value):
    zoomMssg = 'Zoom on reverse primers ±' + str(value) + ' bp.'
    return zoomMssg


@app.callback(
    Output('selected-primers', 'children'),
    [Input('scat', 'clickData'),
     Input('zoomScat', 'clickData')])
def update_output_10(clickForward, clickReverse):
    try:
        frwrdPrimer = clickForward['points'][0]['text']
        # fPosition = clickForward['points'][0]['x']
        # fDegeneracies  = clickForward['points'][0]['y']
        rvrsPrimer = clickReverse['points'][0]['text']
        # rPosition = clickReverse['points'][0]['x']
        # rDegeneracies = clickReverse['points'][0]['y']
        # x = click['points'][0]['x']
        # y = clickData['points'][0]['y']
        # x = str(x)
        # ind = phyloDF[phyloDF.position == x].index.tolist()
        # ind = int(ind[0])
        # print(clickForward)
        # print(clickReverse)
        # Mssg = 'You choose '  + str(frwrdPrimer) + ' as forward primer in ' + str(fPosition) + ' with  ' + (fDegeneracies) + ' and ' + str(rvrsPrimer) + ' as reverse primer in ' + str(rPosition) + ' with  ' + (rDegeneracies)
        Mssg = 'You choose ' + str(frwrdPrimer) + ' --- ' + str(rvrsPrimer)

    except:
        Mssg = ''
    # print(Mssg)
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
        else:
            frwdTemp = frwrdPrimerInfo.TmAp2()
            rvrsTemp = rvrsPrimerInfo.TmAp2()

        tempDist = ff.create_distplot([frwdTemp, rvrsTemp],
                                      ['Forward Primer', 'Reverse Primer'],
                                      curve_type='normal',
                                      show_rug=False,
                                      bin_size=2
                                      )
        tempDist.update_layout(paper_bgcolor='rgba(255,255,255,0.75)'
                               )

    except TypeError:
        tempDist = ff.create_distplot([[0, 0], [0, 0]],
                                      ['Forward Primer', 'Reverse Primer'],
                                      curve_type='normal',
                                      show_rug=False,
                                      bin_size=2
                                      )
        tempDist.update_layout(paper_bgcolor='rgba(255,255,255,0.75)'
                               )

        # print('NO PRIMER CHOOSEN')
    return tempDist


'''
if __name__ == '__main__':
    app.run_server(debug = False)
'''
if __name__ == '__main__':
    app.run_server(debug=False, dev_tools_ui=False, dev_tools_props_check=False)

ind = phyloDF[phyloDF.position == '690-710'].index.tolist()