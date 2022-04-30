import { Box, Stack, TablePagination, Typography } from '@mui/material';
import { AssetTable } from 'asset-list/AssetTable';
import { FieldSpecControl } from 'asset-list/FieldSpecControl';
import { FieldSpec } from 'lib/data-map/view-layers';
import { useCallback, useState } from 'react';
import { ListFeature, useSortedFeatures } from '../asset-list/use-sorted-features';

export const AssetListPage = () => {
  const [page, setPage] = useState(1);
  const [pageSize] = useState(15);

  const layerSpec = {
    layer: 'elec_edges_high', // TODO: enable filtering by sector/subsector/asset type
  };

  const [fieldSpec, setFieldSpec] = useState<FieldSpec>({
    fieldGroup: 'damages_expected',
    fieldDimensions: {
      hazard: 'cyclone',
      rcp: '4.5',
      epoch: 2050,
      protection_standard: 0,
    },
    field: 'ead_mean',
  });

  const { features, pageInfo, loading, error } = useSortedFeatures(layerSpec, fieldSpec, page, pageSize);

  const handleTablePaginationChange = useCallback((event, value) => setPage(value + 1), [setPage]);

  const [selectedFeature, setSelectedFeature] = useState<ListFeature>(null);

  return (
    <>
      <article>
        <Stack spacing={2}>
          <Typography variant="h4">Asset list</Typography>
          <Box>
            <Typography variant="h6">Choose variable to sort assets by</Typography>

            <FieldSpecControl fieldSpec={fieldSpec} onFieldSpec={setFieldSpec} />
          </Box>
          <Stack spacing={2}>
            {loading && <Typography>Loading...</Typography>}
            {error && <Typography>Error: {error.message}</Typography>}
            {!loading && !error && (
              <>
                {/* {count != 0 ? <Pagination count={count} page={page} onChange={handlePaginationChange} /> : null} */}

                <AssetTable
                  features={features}
                  selectedFeature={selectedFeature}
                  onSelectedFeature={setSelectedFeature}
                  footer={
                    pageInfo && (
                      <TablePagination
                        count={pageInfo.total}
                        page={page - 1}
                        onPageChange={handleTablePaginationChange}
                        rowsPerPage={pageSize}
                        rowsPerPageOptions={[pageSize]}
                      />
                    )
                  }
                />
              </>
            )}
          </Stack>
        </Stack>
      </article>
    </>
  );
};
