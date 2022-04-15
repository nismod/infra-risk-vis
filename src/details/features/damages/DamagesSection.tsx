import { Download } from '@mui/icons-material';
import { IconButton, MenuItem, Select, Typography } from '@mui/material';
import { Box } from '@mui/system';
import { HAZARD_DOMAINS } from 'config/hazards/domains';
import { Damage } from 'lib/api-client';
import { downloadFile, titleCase, unique } from 'lib/helpers';
import { useSelect } from 'lib/hooks/use-select';
import _ from 'lodash';
import { useMemo } from 'react';
import { DamageTable } from './DamageTable';
import { EADChart } from './EADChart';

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

function getDamageObject({ hazard, rcp, epoch, mean }: Damage) {
  return {
    key: getDamageKey({ hazard, rcp, epoch }),
    hazard,
    rcp,
    epoch: epoch.toString(),
    ead: mean,
  };
}

function prepareDamages(damages: Damage[]) {
  return damages.filter((d) => d.damage_type === 'direct' && d.protection_standard === 0).map(getDamageObject);
}

function orderDamages(damages: any[]) {
  const lookup = _.fromPairs(damages.map((d) => [d.key, d]));

  return DAMAGES_ORDERING.map(getDamageKey)
    .map((key) => lookup[key])
    .filter(Boolean);
}

function makeDamagesCsv(damages) {
  return 'hazard,rcp,epoch,ead\n' + damages.map((d) => `${d.hazard},${d.rcp},${d.epoch},${d.ead}`).join('\n');
}

export const DamagesSection = ({ fd }) => {
  const damagesData = orderDamages(prepareDamages(fd?.damages ?? []));

  const hazards = useMemo(() => unique(damagesData.map((d) => d.hazard)), [damagesData]);
  const [selectedHazard, setSelectedHazard] = useSelect(hazards);

  const selectedData = useMemo(
    () => (selectedHazard ? damagesData.filter((x) => x.hazard === selectedHazard) : null),
    [selectedHazard, damagesData],
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
        {hazards.length !== 0 && (
          <Select
            variant="standard"
            value={selectedHazard ?? ''}
            onChange={(e) => setSelectedHazard(e.target.value as string)}
          >
            {hazards.map((h) => (
              <MenuItem key={h} value={h}>
                {titleCase(h)}
              </MenuItem>
            ))}
          </Select>
        )}
        {selectedData ? (
          <>
            <Box mt={1}>
              <EADChart
                data={{
                  table: selectedData,
                }}
                actions={false}
                padding={0}
                width={260} // this is currently picked to fit the chart to the sidebar width
                height={150}
                renderer="svg"
              />
            </Box>
            <Box mt={1}>
              <Typography variant="subtitle2">Details</Typography>
              <DamageTable damages={damagesData} />
            </Box>
          </>
        ) : (
          <Typography variant="body2" color="textSecondary">
            No exposure direct damages estimated.
          </Typography>
        )}
      </Box>
    </>
  );
};
