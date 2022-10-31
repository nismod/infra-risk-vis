import { useMemo } from 'react';
import { VegaLite } from 'react-vega';

import { unique } from '@/lib/helpers';

const makeSpec = (yearValues: number[], field_min: string, field: string, field_max: string, field_title: string) => ({
  $schema: 'https://vega.github.io/schema/vega-lite/v5.json',
  data: {
    name: 'table',
  },
  mark: {
    type: 'line',
    point: {
      filled: true,
    },
    tooltip: true,
  },
  encoding: {
    x: {
      field: 'epoch',
      type: 'ordinal',
      title: 'Year',
      axis: {
        gridDash: [2, 2],
        domainColor: '#ccc',
        tickColor: '#ccc',
        values: yearValues,
      },
    },
    y: {
      field: field,
      type: 'quantitative',
      title: field_title,
      axis: {
        gridDash: [2, 2],
        domainColor: '#ccc',
        tickColor: '#ccc',
      },
    },

    color: {
      field: 'rcp',
      type: 'ordinal',
      scale: {
        domain: ['baseline', '2.6', '4.5', '8.5'],
      },
      title: 'RCP',
      legend: {
        orient: 'bottom',
        direction: 'horizontal',
      },
    },
    // the tooltip encoding needs to replicate the field definitions in order to customise their ordering
    tooltip: [
      { field: field, type: 'quantitative', format: ',.3r', title: field_title },
      { field: 'rcp', title: 'RCP' },
      { field: 'epoch', type: 'ordinal', title: 'Year' },
    ],
  },
});

export const ExpectedDamageChart = ({ data, field, field_min, field_max, field_title, ...props }) => {
  const spec = useMemo(
    () => makeSpec(unique<number>(data.table.map((d) => d.epoch)).sort(), field_min, field, field_max, field_title),
    [data, field_min, field, field_max, field_title],
  );

  return <VegaLite data={data} spec={spec as any} {...props} />;
};
