import { useMemo } from 'react';
import { VegaLite } from 'react-vega';

import { unique } from 'lib/helpers';

const makeSpec = (yearValues: number[]) => ({
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
      timeUnit: 'year',
      title: 'Year',
      axis: {
        gridDash: [2, 2],
        domainColor: '#ccc',
        tickColor: '#ccc',
        values: yearValues,
      },
    },
    y: {
      field: 'ead',
      type: 'quantitative',
      title: 'EAD',
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
      { field: 'ead', type: 'quantitative', format: ',.3r', title: 'EAD' },
      { field: 'rcp', title: 'RCP' },
      { field: 'epoch', timeUnit: 'year', title: 'Year' },
    ],
  },
});

export const EADChart = ({ data, ...props }) => {
  const spec = useMemo(() => makeSpec(unique<number>(data.table.map((d) => d.epoch)).sort()), [data]);

  return <VegaLite data={data} spec={spec as any} {...props} />;
};
