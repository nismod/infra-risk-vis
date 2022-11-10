import { Download } from '@mui/icons-material';
import { IconButton } from '@mui/material';
import { Box } from '@mui/system';
import { FC } from 'react';

interface DownloadButtonProps {
  makeContent: () => string;
  filename: string;
  title: string;
  mimeType?: string;
}

export const DownloadButton: FC<DownloadButtonProps> = ({ makeContent, filename, title, mimeType = 'text/csv' }) => {
  return (
    <IconButton title={title} onClick={() => downloadFile(makeContent(), mimeType, filename)}>
      <Download />
    </IconButton>
  );
};

export const ButtonPlacement: FC<{ right?: number }> = ({ children, right = 0 }) => (
  <Box sx={{ position: 'absolute', top: 0, right }}>{children}</Box>
);

// adapted from https://stackoverflow.com/a/44661948/1478817
export function downloadFile(content: string, mimeType: string, fileName: string) {
  const element = document.createElement('a');
  const file = new Blob([content], { type: mimeType });
  element.href = URL.createObjectURL(file);
  element.download = fileName;
  document.body.appendChild(element); // Required for this to work in FireFox
  element.click();
}
