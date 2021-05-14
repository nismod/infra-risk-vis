import React, { Fragment } from 'react';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemText from '@material-ui/core/ListItemText';
import Typography from '@material-ui/core/Typography';
import  { VegaLite } from 'react-vega';

import { titleCase } from './helpers';

const FeatureSidebar = (props) => {
  if (!props.feature) {
    return null
  }
  const f = props.feature.properties;

  return (
    <div className="custom-map-control top-right selected-feature">
      <Typography variant="h6">Selected Asset</Typography>
      <pre style={{display:"none"}}><code>
        {JSON.stringify(f, "", 2)}
      </code></pre>
      <List>
        {
          Object.entries(f).map(([key, value]) => (
            <ListItem>
              <ListItemText
                primary={titleCase(key.replace(/_/g, " "))}
                primaryTypographyProps={{variant:"caption"}}
                secondary={(value ===" ")? "-": value || "-"}
                />
            </ListItem>
          ))
        }
      </List>
      {
        (f.node_id === 'Pump_682')?
          <div style={{paddingLeft:'10px'}}>
          <Typography variant="caption">Stream Flow</Typography>
          <br/>
          <VegaLite
            spec={{
              width: 400,
              height: 200,
              mark: 'line',
              encoding: {
                x: { field: 'date', type: 'temporal', title: 'Date' },
                y: {
                  field: 'flow_1000m3_per_day',
                  type: 'quantitative',
                  title: 'Stream Flow (1000 m³ per day)'
                }
              },
              data: { url: 'flow_cave_river.csv' },
            }}
            // defines the actions available behind the ... menu top-right of chart
            actions={{
              export: true,
              source: false,
              compiled: false,
              editor: false,
            }}
          />
          </div>
        : null
      }
    </div>
  )
}


export default FeatureSidebar;
