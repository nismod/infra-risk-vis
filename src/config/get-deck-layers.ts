import { InteractionTarget } from 'lib/map/interactions/use-interactions';
import { ViewLayer } from 'lib/view-layers';
import { ViewLayerParams } from 'map/get-view-layers-spec';

export function getDeckLayersFunction(
  viewLayers: ViewLayer[],
  viewLayersSpec: Record<string, ViewLayerParams>,
  styleParams: any,
  selectedAsset: InteractionTarget<any>,
) {
  return ({ zoom }: { zoom: number }) =>
    viewLayers.map((layer) => {
      return layer.fn({
        props: { id: layer.id, pickable: !!layer.interactionGroup },
        zoom,
        ...viewLayersSpec[layer.id],
        styleParams,
        selection: selectedAsset?.viewLayer === layer.id ? selectedAsset : undefined,
      });
    });
}
