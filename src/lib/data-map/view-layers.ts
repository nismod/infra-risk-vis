import { DataLoader } from 'lib/data-loader/data-loader';
import { InteractionTarget } from './interactions/use-interactions';

export interface ViewLayerFunctionOptions {
  deckProps: any;
  zoom: number;
  styleParams?: object;
  selection?: InteractionTarget<any>;
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
  spatialType?: string;
  interactionGroup?: string;
  dataManager?: DataManager;
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
