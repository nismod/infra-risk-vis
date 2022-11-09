import reactStringReplace from 'react-string-replace';

export function formatAbbreviations(text, abbreviations: Record<string, string>) {
  for (let [key, value] of Object.entries(abbreviations)) {
    text = reactStringReplace(text, key, (match, i) => (
      <abbr key={match + i} title={value as string}>
        {key}
      </abbr>
    ));
  }

  return text;
}
