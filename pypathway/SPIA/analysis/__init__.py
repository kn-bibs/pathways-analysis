from ..utils import IdMapping
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import math
import pandas as pd
import os
from ..pathviz.utils import plot_json
from copy import deepcopy


class PreProcess:
    @staticmethod
    def deg_stat(data, classes, pos, neg, adjust='fdr_bh'):
        '''
        Basic t-test for certain normalized DataFrame
        If its a RNA SEQ data, use READemption for data process is a better option

        :param data: the pandas dataframe
        :param classes: the class vector
        :param pos: the positive class name
        :param neg: the negative class name
        :param adjust: the multipletest adjust method
        :return: a dataframe contains the result of the basic ttest.
        '''
        data = data.copy()
        PDF = data.groupby(classes, axis=1).get_group(pos)
        CDF = data.groupby(classes, axis=1).get_group(neg)
        ttests = [ttest_ind(PDF.iloc[i], CDF.iloc[i], equal_var=False)[1] for i in range(PDF.shape[0])]
        fc = PDF.mean(axis=1) - CDF.mean(axis=1)
        mul = multipletests(ttests, method=adjust)
        data['fold-change'] = pd.Series(fc, index=data.index)
        data['p-value'] = pd.Series(ttests, index=data.index)
        data['fdr'] = pd.Series(mul[1], index=data.index)
        return data


class EnrichmentResult:
    '''
    This class is the basic class of all enrichment class, contains the necessary result.

    '''

    def __init__(self, df, source_data, target, method, title_index, prop_index, xlabel='-lg(P-Value)'):
        self.df, self.source_data, self.target, self.method, self.xlabel = df, source_data, target, method, xlabel
        self.title_index, self.prop_index = title_index, prop_index
        self.basic_config = {
            'title': {
                'text': 'Enrichment result',
                'subtext': 'Unknown'
            },
            'tooltip': {
                'trigger': 'axis',
                'axisPointer': {
                    'type': 'shadow'
                }
            },
            'grid': {
                'left': '2%',
                'right': '12%',
                'bottom': '3%',
                'containLabel': True
            },
            'xAxis': {
                'type': 'value',
                'boundaryGap': [0, 0.01],
                'name': self.xlabel,
            },
            'yAxis': {
                'type': 'category',
                'data': [],
            },
            'series': [
                {
                    'name': xlabel,
                    'type': 'bar',
                    'data': []
                },
            ]
        }

    def plot(self, count=15, data=False, func=lambda x: -math.log2(x)):
        '''
        Plot the chart in the output area, if it is in the CLI, this function only return an Ipython.display.HTML
        instance

        :param count: the count displayed in the chart, default is 15
        :param data: if True, the config is returned
        :param func: the function to calculate the height of the bar, default is f(x) = -log2(x)
        :return: the Ipython.display.HTML object
        '''
        config = deepcopy(self.basic_config)
        config['yAxis']['data'] = []
        config['series'][0]['data'] = []
        config['title']['subtext'] = self.target
        config['title']['text'] = "{} Enrichment Analysis".format(self.method.upper())
        candidate = []
        for i, x in self.df.iterrows():
            candidate.append([x[self.title_index] if self.title_index >= 0 else i,
                              func(x[self.prop_index] if self.prop_index >= 0 else i)])
        candidate = sorted(candidate, key=lambda x: x[1], reverse=True)
        for x in candidate[:count]:
            config['yAxis']['data'].append({'value': x[0], "textStyle": {"fontSize": 12}})
            config['series'][0]['data'].append(x[1])
        if data:
            return config
        return plot_json(config)

    @property
    def table(self):
        '''
        get the table result of the analysis

        :return: the result DataFrame
        '''
        return self.df

    def overview(self):
        '''
        get the overview of the analysis, contains the method name, arguments and the result stats

        :return: a dict of overview
        '''
        raise NotImplementedError()

    def graph(self):
        '''
        If the result of the enrichment contains the graph relation ship, use this function

        :return: the graph plot
        '''
        raise NotImplementedError()

    def arguments(self):
        raise NotImplementedError()

    def snapshot(self, count=8, func=lambda x: -math.log2(x)):
        '''
        this function return the snapshot for a smaller plot area
        :param count: the count displayed in the chart, default is 15
        :param data: if True, the config is returned
        :param func: the function to calculate the height of the bar, default is f(x) = -log2(x)
        :return: the Ipython.display.HTML object
        '''
        config = deepcopy(self.basic_config)
        config['yAxis']['data'] = []
        config['series'][0]['data'] = []
        config['title']['subtext'] = self.target
        config['title']['text'] = "{} Enrichment Analysis".format(self.method.upper())
        config['color'] = ['#F2F2F2']
        config['xAxis']['lineColor'] = 'transparent'
        config['grid']['top'] = '5%'
        config['xAxis']['splitLine'] = {"show": False}
        config['yAxis']['splitLine'] = {"show": False}
        del config['title']
        candidate = []
        for x in self.df.iterrows():
            candidate.append([x[1][self.title_index], func(x[1][self.prop_index])])
        candidate = sorted(candidate, key=lambda x: x[1], reverse=True)
        for x in candidate[:count]:
            config['yAxis']['data'].append(x[0])
            config['series'][0]['data'].append(x[1])
        return config

    def table_display(self):
        return self.table.to_html()


