import react from '@vitejs/plugin-react';
import fs from 'fs';
import path from 'path';
import { defineConfig } from 'vite';
import svgrPlugin from 'vite-plugin-svgr';
import viteTsconfigPaths from 'vite-tsconfig-paths';

/**
 * To set up a development proxy, create and edit the file dev-proxy/proxy-table.json
 * You can copy and rename one of the included examples.
 * See https://vitejs.dev/config/server-options.html#server-proxy for syntax
 */
let devProxy;
try {
  devProxy = JSON.parse(
    fs.readFileSync(path.resolve(__dirname, 'dev-proxy/proxy-table.json')).toString(),
  );
} catch (err) {
  if (err.code === 'ENOENT') {
    // no dev proxy config, do nothing
  } else {
    console.error('Invalid JSON file: dev-proxy/proxy-table.json');
  }
}

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [react(), viteTsconfigPaths(), svgrPlugin()],
  build: {
    outDir: 'build',
  },
  server: {
    open: true,
    proxy: devProxy,
  },
});
