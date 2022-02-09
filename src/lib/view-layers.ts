import { InteractionTarget } from './map/interactions/use-interactions';

export interface ViewLayerFunctionOptions {
  props: any;
  zoom: number;
  params: any;
  styleParams: object;
  selection: InteractionTarget<any>;
}

export interface ViewLayer {
  id: string;
  fn: (options: ViewLayerFunctionOptions) => any;
  spatialType?: string;
  interactionGroup?: string;
  dataParams?: any;
  getLogicalLayer?: (options: any) => string;
}

export function viewOnlyLayer(id, fn): ViewLayer {
  return {
    id,
    interactionGroup: null,
    fn,
  };
}
