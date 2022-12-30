import { ScaleSequential } from 'd3-scale';
import { ReactNode } from 'react';

import { DataLoader } from '@/lib/data-loader/data-loader';
import { Accessor } from '@/lib/deck/props/getters';

import { InteractionTarget } from './interactions/use-interactions';

export interface FieldSpec {
  fieldGroup: string;
  fieldDimensions?: any;
  field: string;
  fieldParams?: any;
}

export interface ColorSpec {
  scheme: (t: number, n: number) => string;
  scale: (domain: [number, number], interpolator: (t: number, n: number) => string) => ScaleSequential<any, any>;
  range: [number, number];
  empty: string;
  zeroIsEmpty?: boolean;
}
export interface ColorMap {
  fieldSpec: FieldSpec;
  colorSpec: ColorSpec;
}
export interface StyleParams {
  colorMap?: ColorMap;
}
export interface ViewLayerFunctionOptions {
  deckProps: any;
  zoom: number;
  selection?: InteractionTarget<any>;
}

export interface DataManager {
  getDataAccessor: (layer: string, fieldSpec: any) => (d: any) => any;
  getDataLoader: (layer: string, fieldSpec: any) => DataLoader;
}

export interface FormatConfig<D = any> {
  getDataLabel: (fieldSpec: FieldSpec) => string;
  getValueFormatted: (value: D, fieldSpec: FieldSpec) => string | ReactNode;
}

export type ViewLayerDataAccessFunction = (fieldSpec: FieldSpec) => Accessor<any>;
export type ViewLayerDataFormatFunction = (fieldSpec: FieldSpec) => FormatConfig;
export type ViewLayerRenderTooltipFunction = (hover: any) => ReactNode;
export type ViewLayerRenderDetailsFunction = (selection: any) => ReactNode;
export interface ViewLayer<ParamsT = any> {
  id: string;
  params?: ParamsT;
  styleParams?: StyleParams;
  fn: (options: ViewLayerFunctionOptions) => any;
  dataAccessFn?: ViewLayerDataAccessFunction;
  dataFormatsFn?: ViewLayerDataFormatFunction;
  legendDataFormatsFn?: ViewLayerDataFormatFunction;
  spatialType?: string;
  interactionGroup?: string;

  /**
   * new approach for legends: keep the rendering logic inside view layer
   * (currently used for raster layers only)
   */
  renderLegend?: () => ReactNode;
  renderTooltip?: ViewLayerRenderTooltipFunction;
  renderDetails?: ViewLayerRenderDetailsFunction;
}

export function viewOnlyLayer(id, fn): ViewLayer {
  return {
    id,
    interactionGroup: null,
    fn,
  };
}

export interface ViewLayerParams {
  selection?: any;
}
