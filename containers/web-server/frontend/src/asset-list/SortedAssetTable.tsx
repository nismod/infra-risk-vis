import {
  Box,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TablePagination,
  TableRow,
  Typography,
} from '@mui/material';
import React, { FC, ReactNode, useCallback, useEffect, useState } from 'react';

import { FieldSpec } from '@/lib/data-map/view-layers';

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

  const { features, pageInfo, loading, error } = useSortedFeatures(
    layerSpec,
    fieldSpec,
    page,
    pageSize,
  );

  const handleTablePaginationChange = useCallback((event, value) => setPage(value + 1), [setPage]);

  const currentPageFirstItemIndex = (page - 1) * pageSize;

  return (
    <>
      <TableContainer component={Box} height="calc(100% - 48px)" overflow="scroll">
        <Table stickyHeader size="small">
          <TableHead>
            <TableRow>{header}</TableRow>
          </TableHead>

          <TableBody>
            {loading && (
              <TableRow>
                <TableCell colSpan={10} align="center">
                  <Typography variant="body2">Loading...</Typography>
                </TableCell>
              </TableRow>
            )}
            {error && (
              <TableRow>
                <TableCell colSpan={10} align="center">
                  <Typography variant="body2">Error: {error.message}</Typography>
                </TableCell>
              </TableRow>
            )}
            {!loading &&
              !error &&
              features.map((feature, index) =>
                renderRow(feature, index, currentPageFirstItemIndex + index),
              )}
          </TableBody>
        </Table>
      </TableContainer>
      {pageInfo && (
        <TablePagination
          component={Box}
          sx={{
            overflow: 'hidden',
            position: 'absolute',
            bottom: 0,
            width: '100%',
            height: '48px',
          }}
          count={pageInfo.total}
          page={page - 1}
          onPageChange={handleTablePaginationChange}
          rowsPerPage={pageSize}
          rowsPerPageOptions={[pageSize]}
        />
      )}
    </>
  );
};
