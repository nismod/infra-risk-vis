import { FC } from 'react';
import { DefaultValue } from 'recoil';
import { RecoilSync, RecoilSyncOptions } from 'recoil-sync';
import { SetOptional } from 'type-fest';

const DEFAULT_VALUE = new DefaultValue();

const dateRegex =
  // eslint-disable-next-line no-useless-escape
  /^([\+-]?\d{4}(?!\d{2}\b))((-?)((0[1-9]|1[0-2])(\3([12]\d|0[1-9]|3[01]))?|W([0-4]\d|5[0-2])(-?[1-7])?|(00[1-9]|0[1-9]\d|[12]\d{2}|3([0-5]\d|6[1-6])))([T\s]((([01]\d|2[0-3])((:?)[0-5]\d)?|24\:?00)([\.,]\d+(?!:))?)?(\17[0-5]\d([\.,]\d+)?)?([zZ]|([\+-])([01]\d|2[0-3]):?([0-5]\d)?)?)?)?$/;

function parseJSON(text: string) {
  return JSON.parse(text, (key, value) => {
    if (typeof value === 'string' && value.match(dateRegex)) {
      return new Date(value);
    } else return value;
  });
}

const readLocalStorage = (itemKey: string) => {
  const stringVal = localStorage.getItem(itemKey);
  if (stringVal == null) {
    return DEFAULT_VALUE;
  }
  return parseJSON(stringVal);
};

const writeLocalStorage = ({ diff }) => {
  for (const [key, value] of diff) {
    if (value instanceof DefaultValue) {
      localStorage.removeItem(key);
    } else {
      localStorage.setItem(key, JSON.stringify(value));
    }
  }
};

const listenLocalStorage = ({ updateItem, updateAllKnownItems }) => {
  const onStorage = (event: StorageEvent) => {
    // ignore clear() calls
    if (event.storageArea === localStorage && event.key !== null) {
      let parsed: unknown;
      try {
        parsed = event.newValue === null ? DEFAULT_VALUE : parseJSON(event.newValue);
      } catch {
        parsed = DEFAULT_VALUE;
        console.warn({ event }, 'parseJSON failed');
      }

      updateItem(event.key, parsed);
    }
  };

  window.addEventListener('storage', onStorage);

  return () => window.removeEventListener('storage', onStorage);
};

export const RecoilLocalStorageSync: FC<SetOptional<RecoilSyncOptions, 'read' | 'write'>> = (props) => {
  return <RecoilSync read={readLocalStorage} write={writeLocalStorage} listen={listenLocalStorage} {...props} />;
};
