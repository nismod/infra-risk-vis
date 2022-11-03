import { d3Scale, d3ScaleChromatic, discardSides, invertColorScale } from '@/lib/data-map/color-maps';
import { ColorSpec } from '@/lib/data-map/view-layers';

import { DroughtOptionsVariableType, DroughtRiskVariableType } from './metadata';

export const DROUGHT_RISK_COLORMAPS: Record<DroughtRiskVariableType, ColorSpec> = {
  mean_monthly_water_stress_: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateReds,
    range: [0, 15],
    empty: '#ccc',
  },
  epd: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateReds,
    range: [0, 2_500_000],
    empty: '#ccc',
  },
  eael: {
    scale: d3Scale.scaleSequential,
    scheme: invertColorScale(d3ScaleChromatic.interpolateInferno),
    range: [0, 1e9],
    empty: '#ccc',
  },
};

export const DROUGHT_OPTIONS_COLORMAPS: Record<DroughtOptionsVariableType, ColorSpec> = {
  cost_jmd: {
    scale: d3Scale.scaleSequential,
    scheme: discardSides(d3ScaleChromatic.interpolateGreens, 0.2, 0.2),
    range: [0, 1_000_000_000],
    empty: '#ccc',
  },
  population_protected: {
    scale: d3Scale.scaleSequential,
    scheme: discardSides(d3ScaleChromatic.interpolateBlues, 0.2, 0.2),
    range: [0, 100_000],
    empty: '#ccc',
  },
  net_present_value_benefit: {
    scale: d3Scale.scaleSequential,
    scheme: discardSides(d3ScaleChromatic.interpolateBlues, 0.2, 0.2),
    range: [0, 1_000_000_000],
    empty: '#ccc',
  },
  benefit_cost_ratio: {
    scale: d3Scale.scaleSequential,
    scheme: invertColorScale(d3ScaleChromatic.interpolateViridis),
    range: [1, 1.5],
    empty: '#ccc',
  },
};
