"""Write all important results to JSON files
"""
import os
import sys

import pandas as pd
import geopandas as gpd
from collections import defaultdict

def change_depth_string_to_number(x):
    if 'cm' in x:
        return 0.01*float(x.split('cm')[0])
    elif 'm' in x:
        return 1.0*float(x.split('m')[0])
    else:
        return x

def merge_files(mode_edge_dataframe,mode_node_dataframe,value_dataframe,mode_id,value_columns,merge_to_file='edges'):
    if merge_to_file == 'edges':
        mode_dataframe = mode_edge_dataframe
    else:
        mode_dataframe = mode_node_dataframe

    return pd.merge(mode_dataframe,value_dataframe[[mode_id]+value_columns],how='left',on=mode_id)

def flatten_dataframe(df,id_column,anchor_columns,value_columns):
    flat_dict = defaultdict(dict)
    all_cols = []
    for vals in df.itertuples():
        cols = getattr(vals,anchor_columns[0]).replace(' ','_').lower().strip()
        if len(anchor_columns) > 1:
            for a in anchor_columns[1:]: 
                cols += ['_{}'.format(getattr(vals,a).replace(' ','_').lower().strip())][0]
        for v in value_columns:
            flat_dict[getattr(vals,id_column)]['{}_{}'.format(cols,v)] = getattr(vals,v)
            all_cols.append('{}_{}'.format(cols,v))
    
    final_dict = []
    for k,v in flat_dict.items():
        final_dict.append({**{id_column:k},**v})

    del flat_dict
    return pd.DataFrame(final_dict),list(set(all_cols))

def main():
    """Process results
    """
    data_path = os.path.join(os.path.dirname(__file__), '..', 'incoming_data')
    results_path = os.path.join(os.path.dirname(__file__), '..', 'incoming_data','results')
    output_path = os.path.join(os.path.dirname(__file__), '..', 'intermediate_data')

    # Supply input data and parameters
    modes = [
                {
                'sector':'road',
                'edge_file':'road_edges',
                'node_file':'road_nodes',
                'edge_id':'edge_id',
                'node_id':'node_id',
                'edge_attribute_cols':['road_name','road_type','surface',
                                'road_quality','road_service',
                                'min_speed','max_speed',
                                'width','length',
                                'tmda_count'],
                'node_attribute_cols':[],
                'flow_results_cols':['min_total_tons','max_total_tons'],
                'max_flow_colmun':'max_total_tons',
                'fail_results_cols':['min_tr_loss','max_tr_loss',
                                'min_econ_loss','max_econ_loss',
                                'min_econ_impact','max_econ_impact'],
                'flood_results_cols':['min_flood_depth',
                                'max_flood_depth',
                                'min_probability',
                                'max_probability',
                                'min_exposure_length',
                                'max_exposure_length'],
                'risk_results_cols':['ead',
                                'min_eael_per_day',
                                'max_eael_per_day'],
                'adaptation_results_cols':[
                                'options',
                                'ini_adap_cost',
                                'ini_adap_cost_per_km',
                                'tot_maintenance_cost',
                                'tot_maintenance_cost_per_km',
                                'tot_adap_cost',
                                'tot_adap_cost_per_km'],
                },
                {
                'sector':'rail',
                'edge_file':'rail_edges',
                'node_file':'rail_nodes',
                'edge_id':'edge_id',
                'node_id':'node_id',
                'edge_attribute_cols':['operador','linea','length',
                                'min_speed','max_speed'],
                'node_attribute_cols':['nombre','linea','operador'],
                'flow_results_cols':['min_total_tons','max_total_tons'],
                'max_flow_colmun':'max_total_tons',
                'fail_results_cols':['min_tr_loss','max_tr_loss',
                                'min_econ_loss','max_econ_loss',
                                'min_econ_impact','max_econ_impact'],
                'flood_results_cols':['min_flood_depth',
                                'max_flood_depth',
                                'min_probability',
                                'max_probability',
                                'min_exposure_length',
                                'max_exposure_length'],
                'risk_results_cols':[
                                'min_eael_per_day',
                                'max_eael_per_day'],
                'adaptation_results_cols':[],
                },
                {
                'sector':'port',
                'edge_file':'port_edges',
                'node_file':'port_nodes',
                'edge_id':'edge_id',
                'node_id':'node_id',
                'edge_attribute_cols':['length','min_speed','max_speed'],
                'node_attribute_cols':['name','locality','region'],
                'flow_results_cols':['min_total_tons','max_total_tons'],
                'max_flow_colmun':'max_total_tons',
                'fail_results_cols':[],
                'flood_results_cols':['min_flood_depth',
                                'max_flood_depth',
                                'min_probability',
                                'max_probability'],
                'risk_results_cols':[],
                'adaptation_results_cols':[],
                },
                {
                'sector':'air',
                'edge_file':'air_edges',
                'node_file':'air_nodes',
                'edge_id':'edge_id',
                'node_id':'node_id',
                'edge_attribute_cols':['from_iata','to_iata'],
                'node_attribute_cols':['name','iata'],
                'flow_results_cols':['passengers'],
                'max_flow_colmun':'passengers',
                'fail_results_cols':[],
                'flood_results_cols':['min_flood_depth',
                                'max_flood_depth',
                                'min_probability',
                                'max_probability'],
                'risk_results_cols':[],
                'adaptation_results_cols':[],
                },
                {
                'sector':'bridge',
                'edge_file':'bridge_edges',
                'node_file':'bridges',
                'edge_id':'bridge_id',
                'node_id':'bridge_id',
                'edge_attribute_cols':[],
                'node_attribute_cols':[
                                    'pavement_material_asc','substructure_material',
                                    'superstructure_material',
                                    'ruta','structure_type','location',
                                    'width',
                                    'length'],
                'flow_results_cols':['min_total_tons','max_total_tons'],
                'max_flow_colmun':'max_total_tons',
                'fail_results_cols':['min_tr_loss','max_tr_loss',
                                'min_econ_loss','max_econ_loss',
                                'min_econ_impact','max_econ_impact'],
                'flood_results_cols':['min_flood_depth',
                                'max_flood_depth',
                                'min_probability',
                                'max_probability',
                                'min_exposure_length',
                                'max_exposure_length'],
                'risk_results_cols':['ead',
                                'min_eael_per_day',
                                'max_eael_per_day'],
                'adaptation_results_cols':[
                                'options',
                                'ini_adap_cost',
                                'ini_adap_cost_per_km',
                                'tot_maintenance_cost',
                                'tot_maintenance_cost_per_km',
                                'tot_adap_cost',
                                'tot_adap_cost_per_km']
                },

    ]

    geometry_column = 'geometry'
    network_path = os.path.join(data_path,'network')
    for m in modes:
        if m['sector'] == 'road':
            edges_geom = gpd.read_file(os.path.join(network_path,'{}.shp'.format(m['edge_file'])),encoding='utf-8')[[m['edge_id'],geometry_column]]
            edges_csv = pd.read_csv(os.path.join(network_path,'{}.csv'.format(m['edge_file'])),encoding='utf-8-sig')[[m['edge_id']]+m['edge_attribute_cols']]
            edges = pd.merge(edges_geom,edges_csv,how='left',on=[m['edge_id']])
            edges = edges[edges['road_type'].isin(['national','province','rural'])]
            del edges_geom, edges_csv
            nodes = gpd.read_file(os.path.join(network_path,'{}.shp'.format(m['node_file'])),encoding='utf-8')
            merge_type = 'edges'
        if m['sector'] == 'rail':
            edges_geom = gpd.read_file(os.path.join(network_path,'{}.shp'.format(m['edge_file'])),encoding='utf-8')[[m['edge_id'],geometry_column]]
            edges_csv = pd.read_csv(os.path.join(network_path,'{}.csv'.format(m['edge_file'])),encoding='utf-8-sig')[[m['edge_id']]+m['edge_attribute_cols']]
            edges = pd.merge(edges_geom,edges_csv,how='left',on=[m['edge_id']])
            del edges_geom, edges_csv
            nodes = gpd.read_file(os.path.join(network_path,'{}.shp'.format(m['node_file'])),encoding='utf-8')[[m['node_id'],geometry_column]+m['node_attribute_cols']]
            merge_type = 'edges'
        elif m['sector'] in ['air','port']:
            edges = gpd.read_file(os.path.join(network_path,'{}.shp'.format(m['edge_file'])),encoding='utf-8')[[m['edge_id'],geometry_column]+m['edge_attribute_cols']]
            nodes = gpd.read_file(os.path.join(network_path,'{}.shp'.format(m['node_file'])),encoding='utf-8')[[m['node_id'],geometry_column]+m['node_attribute_cols']]
            merge_type = 'edges'
        elif m['sector'] == 'bridge':
            edges = gpd.read_file(os.path.join(network_path,'{}.shp'.format(m['edge_file'])),encoding='utf-8')[[m['edge_id'],geometry_column]]
            nodes_geom = gpd.read_file(os.path.join(network_path,'{}.shp'.format(m['node_file'])),encoding='utf-8')[[m['node_id'],geometry_column]]
            nodes_csv = pd.read_csv(os.path.join(network_path,'{}.csv'.format(m['node_file'])),encoding='utf-8-sig')[[m['node_id']]+m['node_attribute_cols']]
            nodes = pd.merge(nodes_geom,nodes_csv,how='left',on=[m['node_id']])
            del nodes_geom, nodes_csv
            merge_type = 'nodes'

        print ('* Done with reading edges and nodes for {}'.format(m['sector']))

        if merge_type == 'edges':
            mode_df = edges
        else:
            mode_df = nodes
    
        if m['flow_results_cols']:
            if m['sector'] == 'air':
                mode_df = merge_files(mode_df,mode_df,
                        pd.read_csv(os.path.join(data_path,'usage',
                            '{}_passenger.csv'.format(m['sector']))),
                        m['edge_id'],m['flow_results_cols'],merge_to_file=merge_type)
            else:
                mode_df = merge_files(mode_df,mode_df,
                            pd.read_csv(os.path.join(results_path,'flow_mapping_combined',
                                'weighted_flows_{}_100_percent.csv'.format(m['sector']))),
                            m['edge_id'],m['flow_results_cols'],merge_to_file=merge_type)

        print ('* Done with merging flows for {}'.format(m['sector']))

        if m['fail_results_cols']:
            mode_df = merge_files(mode_df,mode_df,
                        pd.read_csv(os.path.join(results_path,
                            'failure_results',
                            'minmax_combined_scenarios',
                            'single_edge_failures_minmax_{}_100_percent_disrupt.csv'.format(m['sector']))),
                        m['edge_id'],m['fail_results_cols'],merge_to_file=merge_type)

        print ('* Done with merging criticality results for {}'.format(m['sector']))

        if m['flood_results_cols'] and m['sector'] in ('road','rail','bridge'):
            f_df,f_cols = flatten_dataframe(pd.read_csv(os.path.join(results_path,
                                        'risk_results',
                                        '{}_hazard_intersections_risk_weights.csv'.format(m['sector']))),
                                        m['edge_id'],
                                        ['hazard_type','climate_scenario'],
                                        m['flood_results_cols'])
            mode_df = merge_files(mode_df,mode_df,
                        f_df,
                        m['edge_id'],
                        f_cols,
                        merge_to_file=merge_type)
            del f_df

        print ('* Done with merging flooding results for {}'.format(m['sector']))

        if m['risk_results_cols']:
            f_df,f_cols = flatten_dataframe(pd.read_csv(os.path.join(results_path,
                                        'risk_results',
                                        '{}_combined_climate_risks.csv'.format(m['sector']))),
                                        m['edge_id'],
                                        ['climate_scenario'],
                                        m['risk_results_cols'])

            mode_df = merge_files(mode_df,mode_df,
                        f_df,
                        m['edge_id'],
                        f_cols,
                        merge_to_file=merge_type)
            del f_df

        print ('* Done with merging risk results for {}'.format(m['sector']))

        if m['adaptation_results_cols']:
            f_df,f_cols = flatten_dataframe(pd.read_csv(os.path.join(results_path,
                                        'adaptation_results',
                                        'combined_climate',
                                        'output_adaptation_{}_costs_fixed_parameters.csv'.format(m['sector']))),
                                        m['edge_id'],
                                        ['climate_scenario'],
                                        m['adaptation_results_cols'])

            mode_df = merge_files(mode_df,mode_df,
                        f_df,
                        m['edge_id'],
                        f_cols,
                        merge_to_file=merge_type)
            del f_df

        print ('* Done with merging adaptation results for {}'.format(m['sector']))

        if merge_type == 'edges':
            edges = mode_df
        else:
            nodes = mode_df

        if m['sector'] in ['air','port']:
            nodes = merge_files(edges,nodes,
                    pd.read_csv(os.path.join(results_path,'network_stats',
                        '{}_ranked_flows.csv'.format(m['sector']))),
                    m['node_id'],m['flow_results_cols'],merge_to_file='nodes')

            nodes = nodes[nodes[m['max_flow_colmun']]>0]
            flood_df = pd.read_csv(os.path.join(results_path,'hazard_scenarios',
                        '{}_hazard_intersections.csv'.format(m['sector'])))
            flood_df['min_depth'] = flood_df.min_depth.apply(
                                                lambda x:change_depth_string_to_number(x))
            flood_df['max_depth'] = flood_df.max_depth.apply(
                                                        lambda x:change_depth_string_to_number(x))
            min_height_prob = flood_df.groupby([m['node_id']] + ['hazard_type','climate_scenario'])['min_depth',
                                                                'probability'].min().reset_index()
            min_height_prob.rename(columns={'min_depth':'min_flood_depth',
                                    'probability': 'min_probability'},inplace=True)
            max_height_prob = flood_df.groupby([m['node_id']] + ['hazard_type','climate_scenario'])['max_depth',
                                                                'probability'].max().reset_index()
            max_height_prob.rename(columns={'max_depth':'max_flood_depth',
                                'probability': 'max_probability'},inplace=True)
            min_max_height_prob = pd.merge(min_height_prob,
                                    max_height_prob,
                                    how='left',
                                    on=[m['node_id'],'hazard_type','climate_scenario'])
            del min_height_prob,max_height_prob

            f_df,f_cols = flatten_dataframe(min_max_height_prob,
                                        m['node_id'],
                                        ['hazard_type','climate_scenario'],
                                        m['flood_results_cols'])
            del min_max_height_prob
            nodes = pd.merge(nodes,f_df,how='left',on=m['node_id'])
            del f_df

        print ('* Done with merging vulnerability results for {}'.format(m['sector']))


        edges.to_file(os.path.join(output_path,'{}.json'.format(m['edge_file'])), driver="GeoJSON", encoding='utf-8')
        nodes.to_file(os.path.join(output_path,'{}.json'.format(m['node_file'])), driver="GeoJSON", encoding='utf-8')

        print ('* Done with writing to JSON for {}'.format(m['sector']))


if __name__ == "__main__":
    main()
