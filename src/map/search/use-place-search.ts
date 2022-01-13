import { useEffect, useState } from 'react';
import { useFetch } from 'use-http';
import { useDebounceCallback } from '@react-hook/debounce';

import { PlaceSearchResult } from './search-state';

interface NominatimSearchResult {
  display_name: string;
  lat: string;
  lon: string;
  boundingbox: string[];
}

function processNominatimData(data: NominatimSearchResult[]): PlaceSearchResult[] {
  return data?.map((x) => ({
    label: x.display_name,
    latitude: parseFloat(x.lat),
    longitude: parseFloat(x.lon),
    boundingBox: {
      minX: parseFloat(x.boundingbox[2]),
      minY: parseFloat(x.boundingbox[0]),
      maxX: parseFloat(x.boundingbox[3]),
      maxY: parseFloat(x.boundingbox[1]),
    },
  }));
}

export function usePlaceSearch(searchValue: string) {
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

  const searchResults: any[] = processNominatimData(data) ?? [];

  return {
    loading,
    error,
    searchResults,
  };
}
