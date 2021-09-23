import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { useCallback, useMemo, useState } from 'react';

import { readPixelsToArray } from '@luma.gl/core';

import { MapParams, useDeckLayers } from './use-map-content';
import { MapViewport } from './MapViewport';
import { MapTooltip } from './tooltip/MapTooltip';
import { FeatureSidebar } from '../FeatureSidebar';
import { FeaturesTooltipContent } from './tooltip/FeaturesTooltipContent';

export const DataMap = ({ background, view, layerSelection }) => {
  const [hoverColor, setHoverColor] = useState<[number, number, number, number]>(null);
  const [hoveredFeatures, setHoveredFeatures] = useState<MapboxGeoJSONFeature[]>([]);
  const [hoverXY, setHoverXY] = useState<[number, number]>(null);

  const [selectedFeatures, setSelectedFeatures] = useState<MapboxGeoJSONFeature[]>([]);

  const mapContentParams = useMemo<MapParams>(
    () => ({
      background,
      view,
      dataLayerSelection: layerSelection,
      highlightedFeature: selectedFeatures?.[0],
    }),
    [background, view, layerSelection, selectedFeatures],
  );

  const onLayerHover = useCallback((info: any) => {
    const { bitmap, object, sourceLayer, x, y } = info;

    if (bitmap || object) {
      setHoverXY([x, y]);
    } else {
      setHoverXY(null);
    }
    if (bitmap) {
      const pixelColor = readPixelsToArray(sourceLayer.props.image, {
        sourceX: bitmap.pixel[0],
        sourceY: bitmap.pixel[1],
        sourceWidth: 1,
        sourceHeight: 1,
        sourceType: undefined,
      });
      if (pixelColor[3]) {
        setHoverColor(pixelColor);
      } else {
        setHoverColor(null);
      }
      setHoveredFeatures(null);
    } else if (object?.properties) {
      setHoverXY([x, y]);
      setHoveredFeatures([object]);
    } else {
      setHoverColor(null);
      setHoveredFeatures(null);
    }
  }, []);

  const onLayerClick = useCallback((info: any) => {
    if (!info.bitmap && info.object) {
      setSelectedFeatures([info.object]);
    } else {
      setSelectedFeatures([]);
    }
  }, []);

  const deckLayers = useDeckLayers(mapContentParams, onLayerHover);

  return (
    <>
      <MapViewport layers={deckLayers} background={background} onHover={onLayerHover} onClick={onLayerClick}>
        <MapTooltip tooltipXY={hoverXY}>
          {hoverColor && (
            <div style={{ backgroundColor: `rgb(${hoverColor[0]},${hoverColor[1]},${hoverColor[2]})` }}>Flooding</div>
          )}
          {hoveredFeatures && <FeaturesTooltipContent features={hoveredFeatures} />}
        </MapTooltip>
      </MapViewport>
      {selectedFeatures.length !== 0 && <FeatureSidebar feature={selectedFeatures[0]} />}
    </>
  );
};
