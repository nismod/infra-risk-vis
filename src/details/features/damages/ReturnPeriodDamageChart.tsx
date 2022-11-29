import { useMemo } from 'react';
import { VegaLite } from 'react-vega';

import { unique } from 'lib/helpers';

const makeSpec = (rpValues: number[], field_min: string, field: string, field_max: string, field_title: string) => ({
  $schema: 'https://vega.github.io/schema/vega-lite/v5.json',
  data: {
    name: 'table',
  },
  layer: [
    {
      mark: 'errorbar',
      encoding: {
        x: {
          field: 'probability',
          type: 'quantitative',
        },
        y: {
          field: field_max,
          title: field_title,
          type: 'quantitative',
        },
        y2: {
          field: field_min,
        },

        color: {
          field: 'rcp',
          type: 'ordinal',
          scale: {
            domain: ['baseline', '4.5', '8.5'],
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
          { field: field, type: 'quantitative', format: ',.3r', title: field_title },
          { field: field_min, type: 'quantitative', format: ',.3r', title: 'Lower bound' },
          { field: field_max, type: 'quantitative', format: ',.3r', title: 'Upper bound' },
          { field: 'rcp', title: 'RCP' },
          { field: 'rp', title: 'Return Period' },
        ],
      },
    },
    {
      mark: {
        type: 'line',
        point: {
          filled: true,
        },
        tooltip: true,
      },
      encoding: {
        x: {
          field: 'probability',
          type: 'quantitative',
          title: 'Probability',
          axis: {
            gridDash: [2, 2],
            domainColor: '#ccc',
            tickColor: '#ccc',
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
            domain: ['baseline', '4.5', '8.5'],
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
          { field: field, type: 'quantitative', format: ',.3r', title: field_title },
          { field: 'rcp', title: 'RCP' },
          { field: 'rp', title: 'Return Period' },
        ],
      },
    },
  ],
});

export const ReturnPeriodDamageChart = ({ data, field, field_min, field_max, field_title, ...props }) => {
  const spec = useMemo(
    () =>
      makeSpec(
        unique<number>(data.table.map((d) => d.rp))
          .sort()
          .reverse(),
        field_min,
        field,
        field_max,
        field_title,
      ),
    [data, field_min, field, field_max, field_title],
  );

  return <VegaLite data={data} spec={spec as any} {...props} />;
};
