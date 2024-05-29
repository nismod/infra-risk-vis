-- Terracotta metadata changes per-database
-- from 0.7.5 to 0.8.3

-- declare version
update metadata set version = '0.8.3';

-- add an integer idx column
alter table key_names add column idx integer;
-- and fill it with numbers from zero in insertion order
with new_key_idx as
    (select k.key_name, (row_number() over () -1) as idx from key_names k)
update key_names k
set idx = (select idx from new_key_idx n where n.key_name = k.key_name);

-- rename filepath column
alter table datasets rename column filepath to path;
