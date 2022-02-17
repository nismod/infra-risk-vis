import _ from 'lodash';

import { LayerDefinition, LayerName, LAYERS } from '../config/layers';
import { ViewName, VIEWS } from '../config/views';

function getSingleViewLayerConfig({ getId, viewLayer }: LayerDefinition) {
  let name: string;
  let dataParams: any;
  if (typeof viewLayer === 'object') {
    dataParams = viewLayer.params;
    name = getId?.(dataParams) ?? viewLayer.baseName;
  } else {
    name = viewLayer;
  }

  return {
    name,
    dataParams,
  };
}

export interface ViewLayerParams {
  visibility: Record<string, boolean>;
  params: object;
  sourceLogicalLayers: string[];
}

export function getViewLayersSpec(
  logicalLayerSelection: Record<string, boolean>,
  view: ViewName,
): Record<string, ViewLayerParams> {
  const viewLayerSpec = {};

  for (const logicalLayerName of VIEWS[view].layers) {
    if (logicalLayerSelection[logicalLayerName] == null) continue;

    const logicalLayerDefinition = LAYERS[logicalLayerName] as LayerDefinition;
    if (logicalLayerDefinition == null) throw new Error(`Logical layer '${logicalLayerName}' is not defined`);

    const { name, dataParams } = getSingleViewLayerConfig(logicalLayerDefinition);

    const singleViewLayerSpec = {
      visibility: {
        [logicalLayerName]: !!logicalLayerSelection[logicalLayerName],
      },
      params: dataParams ?? {},
      sourceLogicalLayers: [logicalLayerName],
    };

    viewLayerSpec[name] = _.merge({}, viewLayerSpec[name], singleViewLayerSpec);
  }

  return viewLayerSpec;
}
