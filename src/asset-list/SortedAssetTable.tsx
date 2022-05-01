import { Box, Table, TableBody, TableContainer, TableHead, TablePagination, TableRow, Typography } from '@mui/material';
import { FieldSpec } from 'lib/data-map/view-layers';
import React, { FC, ReactNode, useCallback, useEffect, useState } from 'react';
import { LayerSpec, ListFeature, useSortedFeatures } from './use-sorted-features';

export const SortedAssetTable: FC<{
  layerSpec: LayerSpec;
  fieldSpec: FieldSpec;
  header: ReactNode;
  renderRow: (feature: ListFeature, localIndex: number, globalIndex: number) => ReactNode;
  pageSize?: number;
}> = ({ layerSpec, fieldSpec, header, renderRow, pageSize = 20 }) => {
  const [page, setPage] = useState(1);

  useEffect(() => {
    setPage(1);
  }, [layerSpec, fieldSpec]);

  const { features, pageInfo, loading, error } = useSortedFeatures(layerSpec, fieldSpec, page, pageSize);

  const handleTablePaginationChange = useCallback((event, value) => setPage(value + 1), [setPage]);

  const currentPageFirstItemIndex = (page - 1) * pageSize;

  return (
    <>
      {loading && <Typography>Loading...</Typography>}
      {error && <Typography>Error: {error.message}</Typography>}
      {!loading && !error && (
        <>
          <TableContainer component={Box}>
            <Table stickyHeader size="small">
              <TableHead>
                <TableRow>{header}</TableRow>
              </TableHead>
              <TableBody sx={{ maxWidth: '100%', maxHeight: '50vh', overflow: 'scroll' }}>
                {features.map((feature, index) => renderRow(feature, index, currentPageFirstItemIndex + index))}
              </TableBody>
            </Table>
          </TableContainer>
          {pageInfo && (
            <TablePagination
              sx={{
                overflow: 'hidden',
              }}
              count={pageInfo.total}
              page={page - 1}
              onPageChange={handleTablePaginationChange}
              rowsPerPage={pageSize}
              rowsPerPageOptions={[pageSize]}
            />
          )}
        </>
      )}
    </>
  );
};
