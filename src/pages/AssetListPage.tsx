import {
  Box,
  FormControl,
  FormLabel,
  List,
  ListItem,
  ListItemText,
  MenuItem,
  Pagination,
  Select,
  Typography,
} from '@mui/material';
import { styled } from '@mui/styles';
import { useMemo, useState } from 'react';
import { useSortedFeatures } from '../asset-list/use-sorted-features';

const ParamDropdown = ({ title, value, onChange, options }) => (
  <FormControl sx={{ minWidth: '10em' }}>
    <FormLabel>{title}</FormLabel>
    <Select value={value} onChange={(e) => onChange(e.target.value)}>
      {options.map((option) => {
        let value, label;
        if (typeof option === 'string' || typeof option === 'number') {
          value = label = option;
        } else {
          value = option.value;
          label = option.label;
        }

        return (
          <MenuItem key={value} value={value}>
            {label}
          </MenuItem>
        );
      })}
    </Select>
  </FormControl>
);
const FieldParamsControl = ({ field, fieldParams, onFieldParams }) => {
  const handleParamChange = (param, value) => {
    onFieldParams({ ...fieldParams, [param]: value });
  };

  return (
    <>
      <Box>
        <ParamDropdown
          title="Hazard"
          value={fieldParams.hazard}
          onChange={(value) => handleParamChange('hazard', value)}
          options={[
            { value: 'total-damages', label: 'All Hazards' },
            { value: 'fluvial', label: 'River Flooding' },
            { value: 'surface', label: 'Surface Flooding' },
            { value: 'coastal', label: 'Coastal Flooding' },
            { value: 'cyclone', label: 'Cyclones' },
          ]}
        />
        <ParamDropdown
          title="RCP"
          value={fieldParams.rcp}
          onChange={(value) => handleParamChange('rcp', value)}
          options={[{ value: 'baseline', label: 'Baseline' }, '2.6', '4.5', '8.5']}
        />
        <ParamDropdown
          title="Epoch"
          value={fieldParams.epoch}
          onChange={(value) => handleParamChange('epoch', value)}
          options={[2010, 2050, 2100]}
        />
      </Box>
    </>
  );
};

const InputSection = styled(Box)({
  marginTop: '1em',
});

interface ListFeature {
  id: number;
  string_id: string;
  bbox_wkt: string;
  value: number;
}

const AssetListItem = ({
  feature,
  onFeatureHover,
}: {
  feature: ListFeature;
  onFeatureHover: (f: ListFeature) => void;
}) => (
  <ListItem
    alignItems="flex-start"
    onMouseOver={(e) => onFeatureHover(feature)}
    onMouseOut={(e) => onFeatureHover(null)}
    sx={{
      '&:hover': {
        backgroundColor: 'rgba(0, 0, 0, 0.05)',
      },
    }}
  >
    <ListItemText
      primary={`Asset ${feature.id} (${feature.string_id})`}
      secondary={
        <>
          <Typography variant="body2">Value: {feature.value}</Typography>
        </>
      }
    />
  </ListItem>
);

export const AssetListPage = () => {
  const [page, setPage] = useState(1);
  const [pageSize] = useState(50);

  const [field, setField] = useState('damages');
  const [fieldParams, setFieldParams] = useState({
    damage_type: 'direct',
    hazard: 'cyclone',
    rcp: '4.5',
    epoch: 2050,
    protection_standard: 0,
  });

  const fieldSpec = useMemo(() => ({ field, fieldParams }), [field, fieldParams]);

  const { features, pageInfo, loading, error } = useSortedFeatures('elec_edges_high', fieldSpec, page, pageSize);

  const count = useMemo(() => (pageInfo ? Math.ceil(pageInfo.total / pageSize) : null), [pageInfo, pageSize]);
  const handlePaginationChange = (event, value) => setPage(value);

  const [hoveredFeature, setHoveredFeature] = useState<ListFeature>(null);

  console.log(hoveredFeature);
  return (
    <>
      <article style={{ position: 'relative', overflow: 'visible' }}>
        <h1>Assets</h1>

        <InputSection>
          <ParamDropdown
            title="Variable"
            value={field}
            onChange={setField}
            options={[{ value: 'damages', label: 'Damages' }]}
          />
        </InputSection>
        <InputSection>
          <Typography variant="subtitle1">Parameters</Typography>
          <FieldParamsControl field={field} fieldParams={fieldParams} onFieldParams={setFieldParams} />
        </InputSection>
        <Box my={2}>
          {loading && <p>Loading...</p>}
          {error && <p>Error: {error.message}</p>}
          {!loading && !error && (
            <>
              {count && <Pagination count={count} page={page} onChange={handlePaginationChange} />}
              <List>
                {features.map((feature) => (
                  <AssetListItem key={feature.id} feature={feature} onFeatureHover={setHoveredFeature} />
                ))}
              </List>
              {count && <Pagination count={count} page={page} onChange={handlePaginationChange} />}
            </>
          )}
        </Box>
        {hoveredFeature && (
          <Box position="sticky" bottom={0} bgcolor="white">
            <Typography variant="h6">Selected asset</Typography>

            <Typography variant="subtitle2">Bounding box</Typography>
            <code>{hoveredFeature.bbox_wkt}</code>
          </Box>
        )}
      </article>
    </>
  );
};
