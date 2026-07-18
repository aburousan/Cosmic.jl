using HTTP
using DataFrames
using CSV

function fetch_gaia_data(query::String)
    url = "https://gea.esac.esa.int/tap-server/tap/sync"
    params = Dict(
        "REQUEST" => "doQuery",
        "LANG" => "ADQL",
        "FORMAT" => "csv",
        "QUERY" => query
    )
    encoded_params = join(["$k=$(HTTP.escapeuri(v))" for (k, v) in params], "&")
    response = HTTP.post(url, ["Content-Type" => "application/x-www-form-urlencoded"], encoded_params)
    if response.status == 200
        data = CSV.read(IOBuffer(response.body), DataFrame)
        return data
    else
        error("Error: ", response.status, "\n", String(response.body))
    end
end

function filter_missing_values(df::DataFrame)
    clean_data = dropmissing(df)
    return clean_data
end

query = """
SELECT TOP 10
    source_id,
    ra,
    dec,
    parallax,
    phot_g_mean_mag
FROM gaiadr3.gaia_source
WHERE ra BETWEEN 0 AND 10
  AND dec BETWEEN -10 AND 10
"""

query1 = """
SELECT
    source_id,
    ra,
    dec,
    parallax,
    phot_g_mean_mag
FROM gaiadr3.gaia_source
WHERE ra BETWEEN 0 AND 10
  AND dec BETWEEN -10 AND 10
"""

function filter_missing_values(df::DataFrame)
    # Filter rows containing missing values
    clean_data = dropmissing(df)
    return clean_data
end


data = fetch_gaia_data(query1)
dat_fil = filter_missing_values(data)
println(dat_fil)