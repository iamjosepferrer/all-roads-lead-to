"""Download Europe's main roads and ferry routes from Overture Maps (OpenStreetMap data) into data_raw/segments.parquet.

Reads the public Overture S3 bucket anonymously with DuckDB; only matching row groups are fetched.
Run: python pipeline/extract.py   (about 5 to 10 minutes, ~490 MB output)
"""
import duckdb, time, sys
from pathlib import Path

RELEASE = "2026-10-21.0"  # Overture release; newer ones are listed at https://docs.overturemaps.org/release/latest/
RAW = Path(__file__).resolve().parent.parent / "data_raw"
RAW.mkdir(exist_ok=True)
OUTFILE = str(RAW / "segments.parquet")
con = duckdb.connect()
con.sql("INSTALL httpfs; LOAD httpfs; INSTALL spatial; LOAD spatial; CREATE SECRET anon (TYPE s3, KEY_ID '', SECRET '', REGION 'us-west-2'); SET threads=4;")
base = f"s3://overturemaps-us-west-2/release/{RELEASE}/theme=transportation/type=segment/*"
t=time.time()
con.sql(f"""COPY (
 SELECT id, subtype, class, subclass, connectors,
        list_transform(speed_limits, lambda s: {{'v': s.max_speed.value, 'u': s.max_speed.unit, 'b': s.between}}) AS spd,
        ST_AsWKB(geometry) AS wkb
 FROM read_parquet('{base}', hive_partitioning=1)
 WHERE bbox.xmin > -25 AND bbox.xmax < 35 AND bbox.ymin > 34 AND bbox.ymax < 72
   AND ((subtype='road' AND class IN ('motorway','trunk','primary')) OR subtype='water')
) TO '{OUTFILE}' (FORMAT parquet, COMPRESSION zstd)""")
print("done", time.time()-t, flush=True)
print(con.sql(f"SELECT subtype, class, count(*) FROM '{OUTFILE}' GROUP BY 1,2").fetchall())
