import { TextField } from '@mui/material';
import { Autocomplete } from '@mui/lab';
import { useState } from 'react';
import { useRecoilState, useSetRecoilState } from 'recoil';

import { usePlaceSearch } from './use-place-search';
import { placeSearchQueryState, placeSearchSelectedResultState } from './search-state';

export const MapSearchField = () => {
  const [searchValue, setSearchValue] = useRecoilState(placeSearchQueryState);
  const [searchResultsOpen, setSearchResultsOpen] = useState(true);

  const { loading, searchResults } = usePlaceSearch(searchValue);

  const setSelectedSearchResult = useSetRecoilState(placeSearchSelectedResultState);

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
        if (reason !== 'createOption') {
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
