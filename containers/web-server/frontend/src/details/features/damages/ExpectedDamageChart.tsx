import _ from 'lodash';
import { FC, useMemo } from 'react';
import { VegaLite } from 'react-vega';

import { unique } from '@/lib/helpers';

import { ExpectedDamageCell } from './ExpectedDamagesSection';

const makeSpec = (yearValues: string[], field_min: string, field: string, field_max: string, field_title: string) => ({
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

// need to map special value to year to maintain chronological ordering on the X axis
function prepareEpoch(epoch: string) {
  return epoch === 'present' ? '2020' : epoch;
}

interface ExpectedDamageChartProps {
  data: {
    table: ExpectedDamageCell[];
  };
  field: keyof ExpectedDamageCell;
  field_min: keyof ExpectedDamageCell;
  field_max: keyof ExpectedDamageCell;
  field_title: string;
}

export const ExpectedDamageChart: FC<ExpectedDamageChartProps> = ({
  data,
  field,
  field_min,
  field_max,
  field_title,
  ...props
}) => {
  // For some reason, Vega was complaining about not being able to extend objects, hence the cloning here.
  // Perhaps it's to do with Recoil freezing state objects
  const clonedData = useMemo(
    () => ({
      table: _.cloneDeep(data.table).map((d) => ({ ...d, epoch: prepareEpoch(d.epoch) })),
    }),
    [data],
  );

  const spec = useMemo(
    () => makeSpec(unique(clonedData.table.map((d) => d.epoch)).sort(), field_min, field, field_max, field_title),
    [clonedData, field_min, field, field_max, field_title],
  );

  return <VegaLite data={clonedData} spec={spec as any} {...props} />;
};
