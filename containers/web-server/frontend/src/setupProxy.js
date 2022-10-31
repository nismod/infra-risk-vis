/**
 * Configure proxies for request to development servers
 *
 * This file must be in node-supported syntax (JS, no ES modules etc.). It does
 * not need to be imported, should be picked up automatically by react-scripts.
 *
 * See guide at https://create-react-app.dev/docs/proxying-api-requests-in-development/#configuring-the-proxy-manually
 */
const { createProxyMiddleware } = require('http-proxy-middleware');

const proxyTable = {
  '/vector': {
    target: 'http://localhost:8080',
    changeOrigin: true,
    pathRewrite: { '^/vector': '/' },
  },

  '/api': {
    target: 'http://localhost:8888',
    changeOrigin: true,
    pathRewrite: { '^/api': '/' },
  },

  // connect to production for frontend-only development

  /*
  '/vector': {
    target: 'https://global.infrastructureresilience.org',
    changeOrigin: true,
    secure: false,
  },

  '/api': {
    target: 'https://global.infrastructureresilience.org',
    changeOrigin: true,
    secure: false,
  },
  */
};

module.exports = function (app) {
  for (const [context, options] of Object.entries(proxyTable)) {
    app.use(context, createProxyMiddleware(options));
  }
};
