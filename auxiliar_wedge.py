import json
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
from pyDatalog import pyDatalog
from pyDatalog.pyDatalog import assert_fact, load, ask

from ipycytoscape import *
import bqplot
from bqplot import pyplot as bqplt
from bqplot import Tooltip

from SPARQLWrapper import SPARQLWrapper, JSON
from math import comb
import os
os.environ["ENDPOINT"]='https://labs.tib.eu/sdm/clarify-kg-7-1/sparql'


def build_query_clarify(input_cui):
    input_cui_uri = ','.join(['<http://clarify2020.eu/entity/' + cui + '>' for cui in input_cui])
    query = """
    select distinct ?EffectorDrugLabel ?AffectedDrugLabel ?Effect ?Impact ?precipitantDrug ?objectDrug ?type
        where {
        {{?s a <http://clarify2020.eu/vocab/DrugDrugInteraction> .  BIND('Pharmacokinetics' as ?type)} 
        UNION {?sim a <http://clarify2020.eu/vocab/SymmetricDrugDrugInteraction> . 
                            ?sim <http://www.w3.org/2000/01/rdf-schema#subClassOf> ?s.BIND('Pharmadynamics' as ?type) }}
        ?s <http://clarify2020.eu/vocab/effect_cui> ?o . 
        ?o <http://clarify2020.eu/vocab/annLabel> ?Effect . 
        ?s <http://clarify2020.eu/vocab/impact> ?Impact .
        ?s <http://clarify2020.eu/vocab/precipitant_drug_cui> ?precipitantDrug .
        ?s <http://clarify2020.eu/vocab/object_drug_cui> ?objectDrug .
        ?precipitantDrug <http://clarify2020.eu/vocab/annLabel> ?EffectorDrugLabel.
        ?objectDrug <http://clarify2020.eu/vocab/annLabel> ?AffectedDrugLabel.

    FILTER (?precipitantDrug in (""" + input_cui_uri + """ ) && ?objectDrug in (""" + input_cui_uri + """))
    }"""
    return query


def store_pharmacokinetic_ddi(effect):
    if effect in ['Excretion_rate', 'Excretory_function', 'Excretion']:
        effect = 'excretion'
    elif effect in ['Process_of_absorption', 'Absorption']:
        effect = 'absorption'
    elif effect in ['Serum_concentration', 'Serum_concentration_of', 'Serum_level', 'Serum_globulin_level', 'Metabolite', 'Active_metabolites']:
        effect = 'serum_concentration'
    elif effect in ['Metabolism']:
        effect = 'metabolism'
    # else:
    #    return 'pharmacodynamic'
    return effect


def rename_impact(impact):
    if impact in ['Increase', 'Higher', 'Worsening']:
        return 'increase'
    return 'decrease'


def get_Labels(input_cui, endpoint):
    labels = {}
    input_cui_uri = ','.join(['<http://clarify2020.eu/entity/' + cui + '>' for cui in input_cui])

    query = """select distinct ?Drug ?drugLabel \n 
    where {?Drug <http://clarify2020.eu/vocab/annLabel> ?drugLabel.\n 
    FILTER (?Drug in (""" + input_cui_uri + """ ))}\n"""

    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    data = results["results"]["bindings"]
    for row in data:
        labels[row["Drug"]["value"].replace("http://clarify2020.eu/entity/", "")] = row["drugLabel"]["value"].lower()

    return list(labels.values())


def query_result_clarify(query, endpoint, labels):
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    dd = {'EffectorDrugLabel': [], 'AffectedDrugLabel': [], 'Effect': [], 'Impact': [], 'precipitantDrug': [],
          'objectDrug': []}
    for r in results['results']['bindings']:
        effect = r['Effect']['value']
        effect = store_pharmacokinetic_ddi(effect)
        dd['Effect'].append(effect.lower())
        impact = r['Impact']['value'].replace('http://clarify2020.eu/entity/', '')
        impact = rename_impact(impact)
        dd['Impact'].append(impact)
        dd['EffectorDrugLabel'].append(r['EffectorDrugLabel']['value'].lower())
        dd['AffectedDrugLabel'].append(r['AffectedDrugLabel']['value'].lower())
        dd['precipitantDrug'].append(r['precipitantDrug']['value'].replace('http://clarify2020.eu/entity/', ''))
        dd['objectDrug'].append(r['objectDrug']['value'].replace('http://clarify2020.eu/entity/', ''))

        if r['Effect']['value']=='Pharmadynamics':
            dd['Effect'].append(effect.lower())
            impact = r['Impact']['value'].replace('http://clarify2020.eu/entity/', '')
            impact = rename_impact(impact)
            dd['Impact'].append(impact)
            dd['EffectorDrugLabel'].append(r['AffectedDrugLabel']['value'].lower())
            dd['AffectedDrugLabel'].append(r['EffectorDrugLabel']['value'].lower())
            dd['precipitantDrug'].append(r['objectDrug']['value'].replace('http://clarify2020.eu/entity/', ''))
            dd['objectDrug'].append(r['precipitantDrug']['value'].replace('http://clarify2020.eu/entity/', ''))

    set_DDIs = pd.DataFrame(dd)
    set_DDIs = set_DDIs.loc[set_DDIs.EffectorDrugLabel.isin(labels)]
    set_DDIs = set_DDIs.loc[set_DDIs.AffectedDrugLabel.isin(labels)]
    set_DDIs.drop_duplicates(keep='first', inplace=True)
    return set_DDIs


def combine_col(corpus, cols):
    # corpus = corpus.apply(lambda x: x.astype(str).str.lower())
    name = '_'.join(cols)
    corpus[name] = corpus[cols].apply(lambda x: '_'.join(x.values.astype(str)), axis=1)
    return corpus


def get_drug_label_by_category(drugs_cui, set_DDIs):
    d_label = set(set_DDIs.loc[set_DDIs.precipitantDrug.isin(drugs_cui)].EffectorDrugLabel.unique())
    d_label.update(set_DDIs.loc[set_DDIs.objectDrug.isin(drugs_cui)].AffectedDrugLabel.unique())
    return d_label


def extract_ddi(onco_drugs, non_onco_drugs, endpoint):
    input_cui = onco_drugs + non_onco_drugs
    labels = get_Labels(input_cui, endpoint)

    query = build_query_clarify(input_cui)
    union = query_result_clarify(query, endpoint, labels)
    union = combine_col(union, ['Effect', 'Impact'])
    set_dsd_label = get_drug_label_by_category(input_cui, union)
    return union, set_dsd_label


def load_data(file):
    onco_drugs = file["Input"]["OncologicalDrugs"]
    non_onco_drugs = file["Input"]["Non_OncologicalDrugs"]
    return extract_ddi(onco_drugs, non_onco_drugs,
                       os.environ["ENDPOINT"])


pyDatalog.create_terms('rdf_star_triple, inferred_rdf_star_triple, wedge, A, B, C, T, T2, wedge_pharmacokinetic')


def build_datalog_model(union):
    pyDatalog.clear()
    for d in union.values:
        # Extensional Database
        assert_fact('rdf_star_triple', d[0], d[1], d[2])
    # Intentional Database
    inferred_rdf_star_triple(A, B, T) <= rdf_star_triple(A, B, T)  # & (T._in(ddiTypeToxicity))
    inferred_rdf_star_triple(A, C, T2) <= inferred_rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (
        T._in(ddiTypeToxicity)) & (T2._in(ddiTypeToxicity)) & (A != C)

    inferred_rdf_star_triple(A, B, T) <= rdf_star_triple(A, B, T)  # & (T._in(ddiTypeEffectiveness))
    inferred_rdf_star_triple(A, C, T2) <= inferred_rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (
        T._in(ddiTypeEffectiveness)) & (T2._in(ddiTypeEffectiveness)) & (A != C)

    wedge(A, B, C, T, T2) <= inferred_rdf_star_triple(A, B, T) & inferred_rdf_star_triple(B, C, T2) & (A != C)

    wedge_pharmacokinetic(A, B, C, T, T2) <= inferred_rdf_star_triple(A, B, T) & inferred_rdf_star_triple(B, C, T2) & (
        T._in(pharmacokinetic_ddi)) & (T2._in(pharmacokinetic_ddi)) & (A != C)


def compute_wedge_datalog(union):
    pyDatalog.clear()
    for d in union.values:
        # Extensional Database
        assert_fact('rdf_star_triple', d[0], d[1], d[2])
    # Intentional Database
    wedge(A, B, C, T, T2) <= rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (A != C)

    wedge_pharmacokinetic(A, B, C, T, T2) <= rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (
        T._in(pharmacokinetic_ddi)) & (T2._in(pharmacokinetic_ddi)) & (A != C)


def dict_wedge_df(dict_wedge, w, drug):
    indirect_ddi = pd.DataFrame(columns=['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact'])
    for i in range(len(w)):
        x = {'EffectorDrugLabel': [w[i][0], drug], 'AffectedDrugLabel': [drug, w[i][1]],
             'Effect_Impact': [w[i][2], w[i][3]]}
        indirect_ddi = pd.concat([indirect_ddi, pd.DataFrame(data=x)])

    indirect_ddi.drop_duplicates(keep='first', inplace=True)
    dict_wedge[drug] = indirect_ddi
    return dict_wedge


def computing_wedge(set_drug_label, ddi_type):
    df_wedge = dict()
    dict_frequency = dict()
    dict_frequency_k = dict()
    max_wedge = len(ddi_type) * comb(len(set_drug_label), 2)
    ddi_k = set.intersection(set(ddi_type), set(pharmacokinetic_ddi))
    max_wedge_k = len(ddi_k) * comb(len(set_drug_label), 2)
    # print(n_ddi, len(set_drug_label), max_wedge)
    for d in set_drug_label:
        w = wedge(A, d, C, T, T2)
        if len(w) > 0:
            df_wedge = dict_wedge_df(df_wedge, w, d)
        dict_frequency[d] = len(w) / max_wedge

        w_k = wedge_pharmacokinetic(A, d, C, T, T2)
        dict_frequency_k[d] = len(w_k) / max_wedge_k
    return df_wedge, dict_frequency, dict_frequency_k


ddiTypeToxicity = ["serum_concentration_increase", "metabolism_decrease", "absorption_increase", "excretion_decrease"]
ddiTypeEffectiveness = ["serum_concentration_decrease", "metabolism_increase", "absorption_decrease",
                        "excretion_increase"]
pharmacokinetic_ddi = ddiTypeToxicity + ddiTypeEffectiveness


def discovering_knowledge(union, set_dsd_label):
    dict_wedge = dict()
    plot_ddi = union[['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact']]
    plot_ddi.drop_duplicates(keep='first', inplace=True)
    build_datalog_model(plot_ddi)
    ddi_type = plot_ddi.Effect_Impact.unique()
    df_wedge, dict_frequency, dict_frequency_k = computing_wedge(set_dsd_label, ddi_type)
    dict_frequency = dict(sorted(dict_frequency.items(), key=lambda item: item[1], reverse=True))
    dict_wedge['DDI_rate'] = dict_frequency
    max_value = max(dict_frequency.values())
    dict_wedge['most_DDI_drug'] = [key for key, value in dict_frequency.items() if value == max_value]

    dict_frequency_k = dict(sorted(dict_frequency_k.items(), key=lambda item: item[1], reverse=True))
    dict_wedge['pharmacokinetic_DDI_rate'] = dict_frequency_k
    max_value = max(dict_frequency_k.values())
    # dict_wedge['most_DDI_drug_pharmacokinetic'] = max(dict_frequency_k, key=dict_frequency_k.get)
    dict_wedge['most_DDI_drug_pharmacokinetic'] = [key for key, value in dict_frequency_k.items() if value == max_value]

    dict_graph_json = dict()
    for k, v in dict_frequency.items():
        if v > 0:
            graph_json = create_json_to_cytoscape(df_wedge[k], k)
            dict_graph_json[k] = graph_json
    return dict_wedge, dict_graph_json


def evaluation_without_deduction(union, set_dsd_label):
    dict_wedge = dict()
    plot_ddi = union[['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact']]
    plot_ddi.drop_duplicates(keep='first', inplace=True)
    compute_wedge_datalog(plot_ddi)
    ddi_type = plot_ddi.Effect_Impact.unique()
    df_wedge, dict_frequency, dict_frequency_k = computing_wedge(set_dsd_label, ddi_type)
    dict_frequency = dict(sorted(dict_frequency.items(), key=lambda item: item[1], reverse=True))
    dict_wedge['DDI_rate'] = dict_frequency
    max_value = max(dict_frequency.values())
    dict_wedge['most_DDI_drug'] = [key for key, value in dict_frequency.items() if value == max_value]

    dict_frequency_k = dict(sorted(dict_frequency_k.items(), key=lambda item: item[1], reverse=True))
    dict_wedge['pharmacokinetic_DDI_rate'] = dict_frequency_k
    max_value = max(dict_frequency_k.values())
    dict_wedge['most_DDI_drug_pharmacokinetic'] = [key for key, value in dict_frequency_k.items() if value == max_value]

    dict_graph_json = dict()
    for k, v in dict_frequency.items():
        if v > 0:
            graph_json = create_json_to_cytoscape(df_wedge[k], k)
            dict_graph_json[k] = graph_json
    return dict_wedge, dict_graph_json


def create_json_to_cytoscape(union, k):
    graph_json = dict()
    graph_json['nodes'] = []
    graph_json['edges'] = []
    drug_id = dict()
    id_x = 0
    for i in range(union.shape[0]):
        precipitant = union.iloc[i]['EffectorDrugLabel']
        object_d = union.iloc[i]['AffectedDrugLabel']
        ddi = union.iloc[i]['Effect_Impact']
        edge = dict()
        edge['data'] = dict()

        if precipitant in drug_id.keys():
            edge['data']['id'] = id_x
            edge['data']['source'] = drug_id[precipitant]
            edge['data']['Effect_Impact'] = ddi
            id_x += 1
        else:
            node = dict()
            node['data'] = dict()
            drug_id[precipitant] = id_x
            node['data']['id'] = id_x
            node['data']['name'] = precipitant
            edge['data']['id'] = id_x + 1
            edge['data']['source'] = id_x
            edge['data']['Effect_Impact'] = ddi
            graph_json['nodes'].append(node)
            id_x += 2
        if object_d in drug_id.keys():
            edge['data']['target'] = drug_id[object_d]
        else:
            node = dict()
            node['data'] = dict()
            drug_id[object_d] = id_x
            node['data']['id'] = id_x
            node['data']['name'] = object_d
            edge['data']['target'] = id_x
            graph_json['nodes'].append(node)
            id_x += 1
            if object_d == k:
                node['classes'] = 'red'  # Single class

        graph_json['edges'].append(edge)

    return graph_json


def create_graph_cytoscape(middle_vertex):
    # load a style dictionary
    with open("styles.json") as fi:
        s = json.load(fi)
    # Create the cytoscape graph widget
    # graph_json = create_json_to_cytoscape(middle_vertex, '')
    cytoscapeobj = CytoscapeWidget()
    cytoscapeobj.graph.add_graph_from_json(middle_vertex, directed=True,
                                           multiple_edges=True)  # , directed=True, input_data['elements']

    cytoscapeobj.set_style(s)
    cytoscapeobj.set_layout(name='breadthfirst', animate=True, nodeSpacing=5)  # concentric,  breadthfirst
    return cytoscapeobj


def get_indirect_ddi(indirect_ddi, dsd):
    deduced_ddi = inferred_rdf_star_triple(C, dsd, T)
    for i in range(len(deduced_ddi)):
        x = {'EffectorDrugLabel': [deduced_ddi[i][0]], 'AffectedDrugLabel': dsd,
             'Effect_Impact': deduced_ddi[i][1]}
        indirect_ddi = pd.concat([indirect_ddi, pd.DataFrame(data=x)])
    return indirect_ddi


# # Whole Graph enriched
def get_graph_enriched(union, set_dsd_label):
    set_DDIs = union[['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact']]
    build_datalog_model(set_DDIs)
    indirect_ddi = pd.DataFrame(columns=['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact'])
    for dsd in set_dsd_label:
        indirect_ddi = get_indirect_ddi(indirect_ddi, dsd)

    graph_ddi = pd.concat([set_DDIs, indirect_ddi])
    graph_ddi.drop_duplicates(keep='first', inplace=True)
    graph_json = create_json_to_cytoscape(graph_ddi, '')
    return graph_ddi, graph_json


def comparision_distribution_wedge(df1, title, x_label, y_label):
    fig = bqplt.figure(title="Wedge frequency distribution", fig_margin={'top':40, 'bottom':40, 'left':60, 'right':0},
                     legend_location="top-right")

    bar_chart  = bqplt.bar(x = list(df1.index), y= [df1["Eval1"], df1["Eval2"]],
                         labels = ["Capturing Knowledge", "Baseline"], display_legend=True)
    bqplt.ylim(0, 1)
    bar_chart.type = "grouped"
    bar_chart.colors = ["tomato", "DeepSkyBlue"]
    bar_chart.tooltip = Tooltip(fields=["x", "y"], labels=[x_label, y_label])
    bqplt.xlabel(x_label)
    bqplt.ylabel(y_label)
    bqplt.show()


if __name__ == '__main__':
    input_list = {
	     "Input":{"OncologicalDrugs":["C3853921"],"Non_OncologicalDrugs":["C0286651","C0012010","C0016410","C0004147"]}
	}
    union, set_dsd_label = load_data(input_list)
    response = discovering_knowledge(union, set_dsd_label)
    r = json.dumps(response, indent=4)
    print(r)