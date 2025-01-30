from pydantic import RootModel


class Polygon(RootModel):
    """Reference to the external GeoJSON Polygon JSON Schema"""

    class Config:
        @staticmethod
        def schema_extra(schema: dict):
            schema.clear()
            schema["$ref"] = "https://geojson.org/schema/Polygon.json"


class MultiPolygon(RootModel):
    """Reference to the external GeoJSON MultiPolygon JSON Schema"""

    class Config:
        @staticmethod
        def schema_extra(schema: dict):
            schema.clear()
            schema["$ref"] = "https://geojson.org/schema/MultiPolygon.json"
