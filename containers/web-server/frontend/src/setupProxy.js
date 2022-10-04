/**
 * Configure proxies for request to development servers
 *
 * This file must be in node-supported syntax (JS, no ES modules etc.). It does
 * not need to be imported, should be picked up automatically by react-scripts.
 *
 * See guide at https://create-react-app.dev/docs/proxying-api-requests-in-development/#configuring-the-proxy-manually
 */
const { createProxyMiddleware } = require('http-proxy-middleware');

module.exports = function (app) {
  app.use(
    '/vector',
    createProxyMiddleware({
      target: 'http://localhost:8080',
      changeOrigin: true,
      pathRewrite: {
        '^/vector': '/', // remove base path
      },
    }),
  );
  app.use(
    '/raster',
    createProxyMiddleware({
      target: 'http://localhost:8888',
      changeOrigin: true,
      pathRewrite: {
        '^/raster': '/', // remove base path
      },
    }),
  );
  app.use(
    '/api',
    createProxyMiddleware({
      target: 'http://localhost:8888',
      changeOrigin: true,
      pathRewrite: {
        '^/api': '/', // remove base path
      },
    }),
  );
};
