import { TextField } from '@material-ui/core';
import { Autocomplete } from '@material-ui/lab';
import { useDebounceCallback } from '@react-hook/debounce';
import { useEffect, useState } from 'react';
import { useRecoilState, useSetRecoilState } from 'recoil';
import { useFetch } from 'use-http';
import { placeSearchQueryState, placeSearchSelectedResultState } from './search-state';

function usePlaceSearch(searchValue: string) {
  const { get, error } = useFetch(
    `https://nominatim.openstreetmap.org/search.php?countrycodes=jm&format=jsonv2&q=${searchValue}`,
  );

  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(false);

  const debouncedGet = useDebounceCallback(async () => {
    try {
      const data = await get();
      setData(data);
    } finally {
      setLoading(false);
    }
  }, 1500);

  useEffect(() => {
    if (searchValue !== '') {
      setLoading(true);
      debouncedGet();
    }
  }, [searchValue, debouncedGet]);

  const searchResults: any[] =
    data?.map((x) => ({ label: x.display_name, latitude: parseFloat(x.lat), longitude: parseFloat(x.lon) })) ?? [];

  return {
    loading,
    error,
    searchResults,
  };
}

export const MapSearchField = () => {
  const [searchValue, setSearchValue] = useRecoilState(placeSearchQueryState);
  const [searchResultsOpen, setSearchResultsOpen] = useState(true);

  const { loading, searchResults } = usePlaceSearch(searchValue);
  // const setSearchResults = useSetRecoilState(placeSearchResultsState);
  const setSelectedSearchResult = useSetRecoilState(placeSearchSelectedResultState);

  // useEffect(() => {
  //   setSearchResults(searchResults);
  // }, [searchResults, setSearchResults]);

  return (
    <Autocomplete
      freeSolo
      openOnFocus
      selectOnFocus
      loading={loading}
      options={loading ? [] : searchResults}
      getOptionLabel={(option) => (typeof option === 'string' ? option : option.label)}
      open={searchResultsOpen}
      onOpen={() => setSearchResultsOpen(true)}
      onClose={() => setSearchResultsOpen(false)}
      onChange={(e, value, reason) => {
        // ignore the change if it's because user pressed enter
        if (reason !== 'create-option') {
          setSelectedSearchResult(value);
        }
      }}
      filterOptions={(x) => x}
      inputValue={searchValue}
      onInputChange={(e, value, reason) => {
        /*
        Ignore the change if it's triggered by the user selecting an option
        This is so that the actual search query stays in the input field
        */
        if (reason !== 'reset') {
          setSearchValue(value);
        }
      }}
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
  );
};
