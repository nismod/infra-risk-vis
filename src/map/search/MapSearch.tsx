import { Box, Button, ButtonGroup, ClickAwayListener, TextField } from '@material-ui/core';
import { Search as SearchIcon } from '@material-ui/icons';
import { Autocomplete } from '@material-ui/lab';
import { useThrottleCallback } from '@react-hook/throttle';
import { useEffect, useState } from 'react';
import useFetch from 'use-http';

function usePlaceSearch(searchValue: string) {
  const { get, error } = useFetch(
    `https://nominatim.openstreetmap.org/search.php?countrycodes=jm&format=jsonv2&q=${searchValue}`,
  );

  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(false);

  const throttledGet = useThrottleCallback(async () => {
    try {
      const data = await get();
      setData(data);
      setLoading(false);
    } catch (err) {
      setLoading(false);
    }
  }, 0.5);

  useEffect(() => {
    if (searchValue !== '') {
      setLoading(true);
      throttledGet();
    }
  }, [searchValue, throttledGet]);

  const searchResults: any[] =
    data?.map((x) => ({ label: x.display_name, latitude: parseFloat(x.lat), longitude: parseFloat(x.lon) })) ?? [];

  return {
    loading,
    error,
    searchResults,
  };
}

export const MapSearch = ({ onPlaceSelect }) => {
  const [expanded, setExpanded] = useState(false);

  const [searchValue, setSearchValue] = useState('');
  const [searchResultsOpen, setSearchResultsOpen] = useState(true);

  const { loading, searchResults } = usePlaceSearch(searchValue);

  return (
    <ClickAwayListener onClickAway={(e) => setExpanded(false)}>
      <Box
        style={{ display: 'flex', flexDirection: 'row' }}
        onKeyDown={(e) => {
          if (e.key === 'Escape') {
            setExpanded(false);
          }
        }}
      >
        <ButtonGroup>
          {/* The ButtonGroup is a hack - causes the button to be narrow. Better find another way */}
          <Button
            variant="contained"
            onClick={() => setExpanded(!expanded)}
            style={{ paddingInline: 0, backgroundColor: 'white' }}
          >
            <SearchIcon />
          </Button>
        </ButtonGroup>
        {expanded && (
          <Autocomplete
            freeSolo
            inputValue={searchValue}
            onInputChange={(e, value) => {
              setSearchValue(value);
            }}
            loading={loading}
            options={searchResults}
            getOptionLabel={(option) => (typeof option === 'string' ? option : option.label)}
            open={searchResultsOpen}
            onOpen={() => setSearchResultsOpen(true)}
            onClose={() => setSearchResultsOpen(false)}
            onChange={(e, value) => {
              if (value) {
                onPlaceSelect({
                  latitude: value.latitude,
                  longitude: value.longitude,
                });
              }
            }}
            filterOptions={(x) => x}
            renderInput={(params) => (
              <TextField
                autoFocus
                variant="standard"
                placeholder="Type place name to search..."
                {...params}
                style={{ width: '300px' }}
              />
            )}
          />
        )}
      </Box>
    </ClickAwayListener>
  );
};
