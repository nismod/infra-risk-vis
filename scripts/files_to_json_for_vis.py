"""Write all important results to JSON files
"""
import os
import sys

import pandas as pd
import geopandas as gpd

def merge_files(mode_edge_file,mode_node_file,value_file_path,mode_id,value_columns,merge_to_file='edges'):
    if merge_to_file == 'edges':
        mode_file = mode_edge_file
    else:
        mode_file = mode_node_file

    value_file = pd.read_csv(value_file_path)[[mode_id]+value_columns]
    return pd.merge(mode_file,value_file,how='left',on=mode_id)


def main():
    """Process results
    """
    data_path = os.path.join(os.path.dirname(__file__), '..', 'incoming_data')
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
                'risk_results_cols':[],
                'vulnerability_results_cols':[],
                'adaptation_results_cols':['hazard_type','climate_scenario', 
                                'min_flood_depth', 'max_flood_depth', 
                                'min_exposure_length', 'max_exposure_length',
                                'min_eael','max_eael',
                                'min_options', 'min_benefit',
                                'min_ini_adap_cost','min_tot_adap_cost','min_bc_ratio',
                                'max_options','max_benefit', 'max_ini_adap_cost','max_tot_adap_cost','max_bc_ratio'],
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
                'vulnerability_results_cols':[],
                'risk_results_cols':['hazard_type','climate_scenario', 
                                'min_flood_depth', 'max_flood_depth', 
                                'min_exposure_length', 'max_exposure_length','risk_wt'],
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
                'vulnerability_results_cols':['hazard_type',
                                    'climate_scenario','probability'],
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
                'flow_results_cols':['passengers_2016'],
                'max_flow_colmun':'passengers_2016',
                'fail_results_cols':[],
                'risk_results_cols':[],
                'vulnerability_results_cols':['hazard_type',
                                    'climate_scenario','probability'],
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
                                    'right_lane_width','left_lane_width','pavement_width_asc','pavement_width_desc',
                                    'length'],
                'flow_results_cols':['min_total_tons','max_total_tons'],
                'max_flow_colmun':'max_total_tons',
                'fail_results_cols':['min_tr_loss','max_tr_loss',
                                'min_econ_loss','max_econ_loss',
                                'min_econ_impact','max_econ_impact'],
                'risk_results_cols':[],
                'vulnerability_results_cols':[],
                'adaptation_results_cols':['hazard_type','climate_scenario', 
                                'min_flood_depth', 'max_flood_depth', 
                                'min_exposure_length', 'max_exposure_length',
                                'min_eael','max_eael',
                                'min_options', 'min_benefit',
                                'min_ini_adap_cost','min_tot_adap_cost','min_bc_ratio',
                                'max_options','max_benefit', 'max_ini_adap_cost','max_tot_adap_cost','max_bc_ratio'],
                },
                
    ]

    geometry_column = 'geometry'
    duration = 10
    network_path = os.path.join(data_path,'network')
    adapt_path = os.path
    for m in modes:
        if m['sector'] == 'road':
            edges_geom = gpd.read_file(os.path.join(network_path,'{}.shp'.format(m['edge_file'])),encoding='utf-8')[[m['edge_id'],geometry_column]]
            edges_csv = pd.read_csv(os.path.join(network_path,'{}.csv'.format(m['edge_file'])),encoding='utf-8-sig')[[m['edge_id']]+m['edge_attribute_cols']]
            edges = pd.merge(edges_geom,edges_csv,how='left',on=[m['edge_id']])
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
            nodes_csv['width'] = nodes_csv['right_lane_width'] + \
                                nodes_csv['left_lane_width'] + \
                                nodes_csv['pavement_width_asc'] + \
                                nodes_csv['pavement_width_desc']
            nodes_csv.drop(['right_lane_width','left_lane_width',
                        'pavement_width_asc','pavement_width_desc'],axis=1,inplace=True)
            nodes = pd.merge(nodes_geom,nodes_csv,how='left',on=[m['node_id']])
            del nodes_geom, nodes_csv
            merge_type = 'nodes'

        print ('* Done with reading edges and nodes for {}'.format(m['sector']))

        if merge_type == 'edges':
            mode_df = edges
        else:
            mode_df = nodes

        if m['flow_results_cols']:
            mode_df = merge_files(mode_df,mode_df,
                        os.path.join(data_path,'flow_results',
                            'weighted_flows_{}_100_percent.csv'.format(m['sector'])),
                        m['edge_id'],m['flow_results_cols'],merge_to_file=merge_type)

        print ('* Done with merging flows for {}'.format(m['sector']))

        if m['fail_results_cols']:
            mode_df = merge_files(mode_df,mode_df,
                        os.path.join(data_path,'failure_results',
                            'single_edge_failures_minmax_{}_100_percent_disrupt.csv'.format(m['sector'])),
                        m['edge_id'],m['fail_results_cols'],merge_to_file=merge_type)

        print ('* Done with merging failure results for {}'.format(m['sector']))

        if m['risk_results_cols']:
            mode_df = merge_files(mode_df,mode_df,
                        os.path.join(data_path,'risk_results',
                            'national_{}_hazard_intersections_risks.csv'.format(m['sector'])),
                        m['edge_id'],m['risk_results_cols'],merge_to_file=merge_type)

            mode_df['min_eael'] = duration*mode_df['risk_wt']*mode_df['min_econ_impact']
            mode_df['max_eael'] = duration*mode_df['risk_wt']*mode_df['max_econ_impact']
            mode_df.drop('risk_wt',axis=1,inplace=True)
            mode_df = mode_df.sort_values(by='max_eael',ascending=False)
            if merge_type == 'edges':
                mode_df.drop_duplicates(subset=[m['edge_id']],keep='first',inplace=True)
            else:
                mode_df.drop_duplicates(subset=[m['node_id']],keep='first',inplace=True)

        print ('* Done with merging risk results for {}'.format(m['sector']))

        if m['adaptation_results_cols']:
            mode_df = merge_files(mode_df,mode_df,
                        os.path.join(data_path,'adaptation_results',
                            '{}_adaptation_summary_{}_days_disruption.csv'.format(m['sector'],duration)),
                        m['edge_id'],m['adaptation_results_cols'],merge_to_file=merge_type)

        print ('* Done with merging adaptation results for {}'.format(m['sector']))

        if merge_type == 'edges':
            edges = mode_df
        else:
            nodes = mode_df

        if m['sector'] in ['air','port']:
            nodes = merge_files(edges,nodes,
                    os.path.join(data_path,'flow_results',
                        '{}_ranked_flows.csv'.format(m['sector'])),
                    m['node_id'],m['flow_results_cols'],merge_to_file='nodes')

            nodes = nodes[nodes[m['max_flow_colmun']]>0]
            flood_df = pd.read_csv(os.path.join(data_path,'vulnerability_results',
                        '{}_vulnerability.csv'.format(m['sector'])))[[m['node_id']]+m['vulnerability_results_cols']]
            flood_df = flood_df.groupby([m['node_id'],'hazard_type','climate_scenario'])['probability'].max().reset_index()
            nodes = pd.merge(nodes,flood_df,how='left',on=m['node_id'])


        print ('* Done with merging vulnerability results for {}'.format(m['sector']))


        edges.to_file(os.path.join(output_path,'{}.json'.format(m['edge_file'])), driver="GeoJSON",encoding='utf-8-sig')
        nodes.to_file(os.path.join(output_path,'{}.json'.format(m['node_file'])), driver="GeoJSON",encoding='utf-8-sig')

        print ('* Done with writing to JSON for {}'.format(m['sector']))


if __name__ == "__main__":
    main()
