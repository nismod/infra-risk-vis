import { Layers, Palette, TableRows } from '@mui/icons-material';
import { SvgIconProps } from '@mui/material';
import { ComponentType } from 'react';

import { MapLegend } from 'map/legend/MapLegend';
import { LayersSidebar } from 'sidebar/SidebarContent';
import { DetailsSidebar } from 'details/DetailsSidebar';

export interface TabConfig {
  id: string;
  label: string;
  IconComponent: ComponentType<SvgIconProps>;
  ContentComponent: ComponentType;
}

/**
 *
 */
export const mobileTabsConfig: TabConfig[] = [
  {
    id: 'layers',
    label: 'Layers',
    IconComponent: Layers,
    ContentComponent: LayersSidebar,
  },
  {
    id: 'legend',
    label: 'Legend',
    IconComponent: Palette,
    ContentComponent: MapLegend,
  },
  {
    id: 'details',
    label: 'Selection',
    IconComponent: TableRows,
    ContentComponent: DetailsSidebar,
  },
];
