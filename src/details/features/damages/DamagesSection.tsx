import { Download } from '@mui/icons-material';
import { FormControl, InputLabel, IconButton, MenuItem, Select, Typography } from '@mui/material';
import { Box } from '@mui/system';
import { HAZARD_DOMAINS } from 'config/hazards/domains';
import { ExpectedDamage } from 'lib/api-client';
import { downloadFile, titleCase, unique } from 'lib/helpers';
import { useSelect } from 'lib/hooks/use-select';
import _ from 'lodash';
import { useMemo } from 'react';
import { DamageTable } from './DamageTable';
import { ExpectedDamageChart } from './ExpectedDamageChart';

const DAMAGES_ORDERING = (() => {
  const ordering = [];
  for (const [hazard, hazardDomain] of Object.entries(HAZARD_DOMAINS)) {
    for (const rcp of hazardDomain.paramDomains.rcp) {
      for (const epoch of hazardDomain.paramDomains.epoch) {
        ordering.push({
          hazard,
          rcp,
          epoch,
        });
      }
    }
  }
  return ordering;
})();

function getDamageKey({ hazard, rcp, epoch }) {
  return `${hazard}__rcp_${rcp}__epoch_${epoch}__conf_None`;
}

interface ExpectedDamageCell {
  key: string;
  hazard: string;
  rcp: string;
  epoch: string;
  ead_mean: number,
  ead_amin: number,
  ead_amax: number,
  eael_mean: number,
  eael_amin: number,
  eael_amax: number,
}
function getExpectedDamageObject({
  hazard,
  rcp,
  epoch,
  ead_mean,
  ead_amin,
  ead_amax,
  eael_mean,
  eael_amin,
  eael_amax
}: ExpectedDamage): ExpectedDamageCell {
  return {
    key: getDamageKey({ hazard, rcp, epoch }),
    hazard,
    rcp,
    epoch: epoch.toString(),
    ead_mean,
    ead_amin,
    ead_amax,
    eael_mean,
    eael_amin,
    eael_amax,
  };
}

function prepareExpectedDamages(expectedDamages: ExpectedDamage[]) {
  return expectedDamages.filter((d) => d.protection_standard === 0).map(getExpectedDamageObject);
}

function orderDamages(damages: ExpectedDamageCell[]) {
  const lookup = _.fromPairs(damages.map((d) => [d.key, d]));

  return DAMAGES_ORDERING.map(getDamageKey)
    .map((key) => lookup[key])
    .filter(Boolean);
}

function makeDamagesCsv(damages: ExpectedDamageCell[]) {
  return 'hazard,rcp,epoch,ead_mean,ead_amin,ead_amax,eael_mean,eael_amin,eael_amax\n' + damages.map(
    (d) => `${d.hazard},${d.rcp},${d.epoch},${d.ead_mean},${d.ead_amin},${d.ead_amax},${d.eael_mean},${d.eael_amin},${d.eael_amax}`).join('\n');
}

export const DamagesSection = ({ fd }) => {
  const damagesData = orderDamages(prepareExpectedDamages(fd?.damages_expected ?? []));

  const hazards = useMemo(() => unique(damagesData.map((d) => d.hazard)), [damagesData]);
  const [selectedHazard, setSelectedHazard] = useSelect(hazards);

  const selectedData = useMemo(
    () => (selectedHazard ? damagesData.filter((x) => x.hazard === selectedHazard) : null),
    [selectedHazard, damagesData],
  );

  const has_eael = useMemo(
    () => (selectedData ? selectedData.some((d)=> d.eael_amax > 0) : null),
    [selectedData],
  );

  return (
    <>
      <Box py={2}>
        <Box position="relative">
          <Typography variant="h6">Risk</Typography>
          <IconButton
            sx={{
              position: 'absolute',
              top: 0,
              right: 0,
            }}
            title="Download CSV with damages data"
            onClick={() => downloadFile(makeDamagesCsv(damagesData), 'text/csv', `feature_${fd.id}_damages.csv`)}
          >
            <Download />
          </IconButton>
        </Box>
        {hazards.length ? (
          hazards.length === 1?
            (
              <FormControl fullWidth sx={{my:2}} disabled>
                <InputLabel>Hazard</InputLabel>
                <Select
                  value={selectedHazard ?? ''}
                  onChange={(e) => setSelectedHazard(e.target.value as string)}
                  >
                  {hazards.map((h) => (
                    <MenuItem key={h} value={h}>
                      {titleCase(h)}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
            ) : (
              <FormControl fullWidth sx={{my:2}}>
                <InputLabel>Hazard</InputLabel>
                <Select
                  value={selectedHazard ?? ''}
                  onChange={(e) => setSelectedHazard(e.target.value as string)}
                  >
                  {hazards.map((h) => (
                    <MenuItem key={h} value={h}>
                      {titleCase(h)}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
            )
          ) : null
        }
        {selectedData ? (
          <>
            <Box mt={1}>
              <Typography variant="subtitle2">Expected Annual Damages</Typography>
              <ExpectedDamageChart
                data={{
                  table: selectedData,
                }}
                field='ead_mean'
                field_min='ead_amin'
                field_max='ead_amax'
                field_title='EAD (J$)'
                actions={false}
                padding={0}
                width={260} // this is currently picked to fit the chart to the sidebar width
                height={150}
                renderer="svg"
              />
            </Box>
            {has_eael? (
              <Box mt={1}>
                <Typography variant="subtitle2">Expected Annual Economic Losses</Typography>
                <ExpectedDamageChart
                  data={{
                    table: selectedData,
                  }}
                  field='eael_mean'
                  field_min='eael_amin'
                  field_max='eael_amax'
                  field_title='EAEL (J$/day)'
                  actions={false}
                  padding={0}
                  width={260} // this is currently picked to fit the chart to the sidebar width
                  height={150}
                  renderer="svg"
                  />
              </Box>
              ) : null
            }
            <Box mt={1}>
              <DamageTable damages={selectedData} />
            </Box>
          </>
        ) : (
          <Typography variant="body2" color="textSecondary">
            No direct damages or indirect losses estimated.
          </Typography>
        )}
      </Box>
    </>
  );
};
