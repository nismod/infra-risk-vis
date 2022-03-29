import { InteractionTarget } from './interactions/use-interactions';

export interface ViewLayerFunctionOptions {
  deckProps: any;
  zoom: number;
  styleParams?: object;
  selection?: InteractionTarget<any>;
}

export interface ViewLayer {
  id: string;
  params?: any;
  group: string;
  fn: (options: ViewLayerFunctionOptions) => any;
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
