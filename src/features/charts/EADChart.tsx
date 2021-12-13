import { createClassFromSpec } from 'react-vega';

export const EADChart = createClassFromSpec({
  spec: {
    $schema: 'https://vega.github.io/schema/vega-lite/v5.json',
    data: {
      name: 'table',
    },
    mark: {
      type: 'line',
      point: {
        filled: true,
      },
    },
    encoding: {
      x: { field: 'epoch', timeUnit: 'year', title: 'Year' },
      y: { field: 'ead', type: 'quantitative', title: 'EAD' },
      color: {
        field: 'rcp',
        type: 'nominal',
        title: 'RCP',
        legend: {
          orient: 'bottom',
          direction: 'horizontal',
        },
      },
    },
  },
});
