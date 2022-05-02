import { useMemo } from 'react';
import { VegaLite } from 'react-vega';

import { unique } from 'lib/helpers';

const makeSpec = (yearValues: number[], field_key: string, field_title: string) => ({
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
      field: field_key,
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
        // Could do custom colours
        // range: ["#e7ba52", "#c7c7c7", "#aec7e8", "#1f77b4"]
      },
      title: 'RCP',
      legend: {
        orient: 'bottom',
        direction: 'horizontal',
      },
    },
    // the tooltip encoding needs to replicate the field definitions in order to customise their ordering
    tooltip: [
      { field: field_key, type: 'quantitative', format: ',.3r', title: field_title },
      { field: 'rcp', title: 'RCP' },
      { field: 'epoch', timeUnit: 'year', title: 'Year' },
    ],
  },
});

export const ExpectedDamageChart = ({ data, field_key, field_title, ...props }) => {
  const spec = useMemo(
    () => makeSpec(
      unique<number>(data.table.map((d) => d.epoch)).sort(),
      field_key,
      field_title
    ),
    [data]
  );

  return <VegaLite data={data} spec={spec as any} {...props} />;
};
