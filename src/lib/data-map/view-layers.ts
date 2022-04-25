import { DataLoader } from 'lib/data-loader/data-loader';
import { InteractionTarget } from './interactions/use-interactions';

export interface FieldSpec {
  fieldGroup: string;
  fieldDimensions?: any;
  field: string;
}

export interface ColorMap {
  colorField: FieldSpec;
  colorScheme: string;
}
export interface StyleParams {
  colorMap?: ColorMap;
}
export interface ViewLayerFunctionOptions {
  deckProps: any;
  zoom: number;
  styleParams?: StyleParams;
  selection?: InteractionTarget<any>;
}

export interface ViewLayerDataFunctionOptions {
  styleParams?: StyleParams;
}

export interface DataManager {
  getDataAccessor: (layer: string, fieldSpec: any) => (d: any) => any;
  getDataLoader: (layer: string, fieldSpec: any) => DataLoader;
}
export interface DataAccess {
  dataAccessor: (d: any) => any;
  dataLoader: DataLoader;
}
export interface ViewLayer {
  id: string;
  params?: any;
  group: string;
  fn: (options: ViewLayerFunctionOptions) => any;
  dataAccessFn?: (options: ViewLayerDataFunctionOptions) => DataAccess;
  spatialType?: string;
  interactionGroup?: string;
}

export function viewOnlyLayer(id, fn): ViewLayer {
  return {
    id,
    group: null,
    interactionGroup: null,
    fn,
  };
}

export interface ViewLayerParams {
  selection?: any;
  styleParams?: any;
}
