import { Download } from '@mui/icons-material';
import { FormControl, InputLabel, IconButton, MenuItem, Select, Typography } from '@mui/material';
import { Box } from '@mui/system';
import { HAZARD_DOMAINS } from 'config/hazards/domains';
import { ExpectedDamage, ReturnPeriodDamage } from 'lib/api-client';
import { downloadFile, titleCase, unique } from 'lib/helpers';
import { useSelect } from 'lib/hooks/use-select';
import _ from 'lodash';
import { useMemo } from 'react';
import { DamageTable } from './DamageTable';
import { RPDamageTable } from './RPDamageTable';
import { ExpectedDamageChart } from './ExpectedDamageChart';
import { ReturnPeriodDamageChart } from './ReturnPeriodDamageChart';

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

const RP_ORDERING = (() => {
  const ordering = [];
  for (const [hazard, hazardDomain] of Object.entries(HAZARD_DOMAINS)) {
    for (const rp of hazardDomain.paramDomains.returnPeriod) {
      for (const rcp of hazardDomain.paramDomains.rcp) {
        for (const epoch of hazardDomain.paramDomains.epoch) {
          ordering.push({
            hazard,
            rcp,
            epoch,
            rp
          });
        }
      }
    }
  }
  return ordering;
})();

function getDamageKey({ hazard, rcp, epoch }) {
  return `${hazard}__rcp_${rcp}__epoch_${epoch}__conf_None`;
}
function getRPDamageKey({ hazard, rcp, epoch, rp }) {
  return `${hazard}__rcp_${rcp}__epoch_${epoch}__rp_${rp}__conf_None`;
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
function getExpectedDamageObject(d: ExpectedDamage): ExpectedDamageCell {
  return {
    key: getDamageKey({ hazard: d.hazard, rcp: d.rcp, epoch: d.epoch }),
    hazard: d.hazard,
    rcp: d.rcp,
    epoch: d.epoch.toString(),
    ead_mean: d.ead_mean,
    ead_amin: d.ead_amin,
    ead_amax: d.ead_amax,
    eael_mean: d.eael_mean,
    eael_amin: d.eael_amin,
    eael_amax: d.eael_amax,
  };
}

interface RPDamageCell {
  key: string;
  hazard: string;
  rcp: string;
  rp: number;
  probability: number;
  epoch: string;
  damage_mean: number,
  damage_amin: number,
  damage_amax: number,
  loss_mean: number,
  loss_amin: number,
  loss_amax: number,
}

function getRPDamageObject(d: ReturnPeriodDamage): RPDamageCell {
  return {
    key: getRPDamageKey({ hazard: d.hazard, rcp: d.rcp, epoch: d.epoch, rp: d.rp }),
    hazard: d.hazard,
    rcp: d.rcp,
    rp: d.rp,
    probability: 1 / d.rp,
    epoch: d.epoch.toString(),
    damage_mean: d.damage_mean,
    damage_amin: d.damage_amin,
    damage_amax: d.damage_amax,
    loss_mean: d.loss_mean,
    loss_amin: d.loss_amin,
    loss_amax: d.loss_amax,
  };
}

function prepareExpectedDamages(expectedDamages: ExpectedDamage[]) {
  return expectedDamages.filter((d) => d.protection_standard === 0).map(getExpectedDamageObject);
}

function prepareRPDamages(rpDamages: ReturnPeriodDamage[]) {
  return rpDamages.map(getRPDamageObject);
}

function orderDamages(damages: ExpectedDamageCell[]) {
  const lookup = _.fromPairs(damages.map((d) => [d.key, d]));

  return DAMAGES_ORDERING.map(getDamageKey)
    .map((key) => lookup[key])
    .filter(Boolean);
}

function orderRPDamages(damages: RPDamageCell[]) {
  const lookup = _.fromPairs(damages.map((d) => [d.key, d]));

  return RP_ORDERING.map(getRPDamageKey)
    .map((key) => lookup[key])
    .filter(Boolean);
}

function makeDamagesCsv(damages: ExpectedDamageCell[]) {
  return 'hazard,rcp,epoch,ead_mean,ead_amin,ead_amax,eael_mean,eael_amin,eael_amax\n' + damages.map(
    (d) => `${d.hazard},${d.rcp},${d.epoch},${d.ead_mean},${d.ead_amin},${d.ead_amax},${d.eael_mean},${d.eael_amin},${d.eael_amax}`).join('\n');
}

function makeRPDamagesCsv(damages: RPDamageCell[]) {
  return 'hazard,rcp,epoch,rp,damage_mean,damage_amin,damage_amax,loss_mean,loss_amin,loss_amax\n' + damages.map(
    (d) => `${d.hazard},${d.rcp},${d.epoch},${d.rp},${d.damage_mean},${d.damage_amin},${d.damage_amax},${d.loss_mean},${d.loss_amin},${d.loss_amax}`).join('\n');
}

export const DamagesSection = ({ fd }) => {
  const damagesData = orderDamages(prepareExpectedDamages(fd?.damages_expected ?? []));
  const returnPeriodDamagesData = orderRPDamages(prepareRPDamages(fd?.damages_return_period ?? []));

  const hazards = useMemo(() => unique(damagesData.map((d) => d.hazard)), [damagesData]);
  const epochs = useMemo(() => unique(damagesData.map((d) => d.epoch)).sort(), [damagesData]);
  const [selectedHazard, setSelectedHazard] = useSelect(hazards);
  const [selectedEpoch, setSelectedEpoch] = useSelect(epochs);

  const selectedData = useMemo(
    () => (selectedHazard ? damagesData.filter((x) => x.hazard === selectedHazard) : null),
    [selectedHazard, damagesData],
  );
  const selectedRPData = useMemo(
    () => (selectedHazard ? returnPeriodDamagesData.filter((x) => x.hazard === selectedHazard && x.epoch === selectedEpoch) : null),
    [selectedHazard, selectedEpoch, returnPeriodDamagesData],
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
          <FormControl fullWidth sx={{my:2}} disabled={hazards.length === 1}>
            <InputLabel>Hazard</InputLabel>
            <Select
              label="Hazard"
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
      <Box py={2}>
        <Box position="relative">
          <Typography variant="h6">Return Period Damages</Typography>
          <IconButton
            sx={{
              position: 'absolute',
              top: 0,
              right: 0,
            }}
            title="Download CSV with return period data"
            onClick={() => downloadFile(makeRPDamagesCsv(returnPeriodDamagesData), 'text/csv', `feature_${fd.id}_damages_rp.csv`)}
          >
            <Download />
          </IconButton>
        </Box>
        {epochs.length ? (
          <FormControl fullWidth sx={{my:2}} disabled={epochs.length === 1}>
            <InputLabel>Epoch</InputLabel>
            <Select
              label="Epoch"
              value={selectedEpoch ?? ''}
              onChange={(e) => setSelectedEpoch(e.target.value as string)}
              >
              {epochs.map((h) => (
                <MenuItem key={h} value={h}>
                  {titleCase(h)}
                </MenuItem>
              ))}
            </Select>
          </FormControl>
          ) : null
        }
        {selectedRPData ? (
          <>
            <Box mt={1}>
              <ReturnPeriodDamageChart
                data={{
                  table: selectedRPData,
                }}
                field_key='damage_mean'
                field_title='Damage (J$)'
                actions={false}
                padding={0}
                width={260} // this is currently picked to fit the chart to the sidebar width
                height={150}
                renderer="svg"
              />
            </Box>
            <Box mt={1}>
              <RPDamageTable damages={selectedRPData} />
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
