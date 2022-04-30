import { TablePagination, Typography } from '@mui/material';
import { FieldSpec } from 'lib/data-map/view-layers';
import { FC, useCallback, useState } from 'react';
import { AssetTable } from './AssetTable';
import { LayerSpec, ListFeature, useSortedFeatures } from './use-sorted-features';

export const SortedAssetTable: FC<{
  layerSpec: LayerSpec;
  fieldSpec: FieldSpec;
  selectedFeature?: ListFeature;
  onSelectedFeature?: (feature: ListFeature) => void;
  hoveredFeature?: ListFeature;
  onHoveredFeature?: (feature: ListFeature) => void;
}> = ({ layerSpec, fieldSpec, selectedFeature, onSelectedFeature, hoveredFeature, onHoveredFeature }) => {
  const [page, setPage] = useState(1);
  const [pageSize] = useState(15);

  const { features, pageInfo, loading, error } = useSortedFeatures(layerSpec, fieldSpec, page, pageSize);

  const handleTablePaginationChange = useCallback((event, value) => setPage(value + 1), [setPage]);

  return (
    <>
      {loading && <Typography>Loading...</Typography>}
      {error && <Typography>Error: {error.message}</Typography>}
      {!loading && !error && (
        <>
          <AssetTable
            features={features}
            selectedFeature={selectedFeature}
            onSelectedFeature={onSelectedFeature}
            hoveredFeature={hoveredFeature}
            onHoveredFeature={onHoveredFeature}
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
    </>
  );
};
