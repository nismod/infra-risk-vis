const _ = require('lodash');
const TsconfigPathsPlugin = require('tsconfig-paths-webpack-plugin');

module.exports = function override(config, env) {
  config = _.merge(config, {
    resolve: {
      plugins: [new TsconfigPathsPlugin()],
    },
  });

  return config;
};
