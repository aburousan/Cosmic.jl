### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ eb918f38-b7cc-11ee-1849-e9a292ae285b
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

    using Cosmic
end

# ╔═╡ 1473e933-d04a-4a0b-90d8-ab33df9ff3ee
using Plots,LaTeXStrings

# ╔═╡ 8d34d1a3-6a28-461c-9e03-cc047654ce84
md"Now we can use Cosmic Package just use."

# ╔═╡ ebf05ec9-dc2e-4ac9-a231-fd05e170e0a2
md"The first step is to create your Cosmology system. I will be using Flat model."

# ╔═╡ 667c6bda-f78c-4923-b93c-1942f5b1f258
c = cosmology()#by default values will be taken.

# ╔═╡ dc415649-5683-415e-a1b8-a31f48080d94
md"The first value is h, then comes $\Omega_\Lambda$, then $\Omega_m$ and finally $\Omega_r$."

# ╔═╡ 0314d0ad-d666-41fb-8e5a-bd4a0a497ba6
cop = cosmology(Ωk =0.1)#For some other model can be used

# ╔═╡ a2a0d4d9-46cc-4f18-9311-de22ea571b4e
md"For conversion between redshift and scale factor, there are two functions."

# ╔═╡ 361cb06f-4e96-48a1-a23e-e38106ef4ed8
a(100)#computes the scalefactor for z=100

# ╔═╡ 52930937-1f73-470d-b041-9558a11a7658
redshift(0.0099099)#compute the redshift for a = 0.0099099

# ╔═╡ 5d11e724-326a-436f-94f4-211f99d30990
md"To get the hubble's constant at some redshift use H(c::AbstractCosmology,z) as"

# ╔═╡ 7e4c3515-02f7-498b-829f-23e78997fb70
H(c,100)#takes your cosmology model and z

# ╔═╡ 3aa50f1f-04e8-4d0d-a650-c0ef8e75089f
md"Let's try plotting $\frac{H(z)}{(1+z)}$ vs $z$."

# ╔═╡ 2ee6b0bf-3511-4fea-8132-9896f5f112ec
begin
	z_vals = range(0, 2.5, length=1000)
	H_1z(z) = H(c,z)/(1+z)
	H_vals = H_1z.(z_vals)
	plot(z_vals,H_vals,lw=2.5, label="")
	xlabel!(L"Redshift ($z$)")
	ylabel!(L"$\frac{H(z)}{1+z}$")
end

# ╔═╡ 4984ac71-b8eb-4085-a20a-9a2ac5f439a4
md"We can also compute Hubble's distance. For that use hubble_distance(c,z) function."

# ╔═╡ 2bb09fc2-f284-4d1e-8d58-55365d31dc20
begin
	# z_vals = range(0, 2.5, length=1000)
	Hd_1z(z) =  hubble_distance(c,z)
	Hd_vals = Hd_1z.(z_vals)
	plot(z_vals,Hd_vals,lw=2.5, label="")
	xlabel!(L"Redshift ($z$)")
	ylabel!(L"$H_d(z)$")
end

# ╔═╡ d3e187e5-13ee-4edb-bfa0-25aee07b712e
md"We can also find the age of the universe at some certain z value."

# ╔═╡ b02a5df0-b196-4418-84b0-52633d6e49e0
age(c,0)#z = 0 means today's time.

# ╔═╡ 511ade9b-ccf8-4edd-a5c5-edf2ba869815
md"This gives us the age of the universe in sec. If you want to get the age in Gigayear, just convert it or else use ageGyr(c,z)."

# ╔═╡ 7ed40ab7-c190-4bf0-bcf6-eccbb422784e
ageGyr(c,0)# age of our universe.

# ╔═╡ ac2675fd-b0c9-4ba7-95eb-79f72bddc678
md"For computing the scale factor as a function of time, there are two functions, scalefact_part and scalefact. The first one is efficient but it second one is more accurate. Both takes input in Gyr."

# ╔═╡ f1b373f4-43bd-4312-bcaf-49ff3fc2d31a
scalefact_part(c,100)

# ╔═╡ efabfaa1-f0bc-466a-bb98-49795c44e200
scale_fact(c,100)

# ╔═╡ b9dc315c-5a35-41fd-b27d-f1f9a4fc7e81
md"Let's make a plot for both of them."

# ╔═╡ 1e09ba10-91d6-46f0-bc8f-b9b1a82220fd
begin
	tgyr = range(0.1,30.8,length=10_000)
	dgf(t) = scale_fact(c,t)
	a_vals = dgf.(tgyr)
	plot(tgyr,a_vals,lw=2.5, label="using scale_fact")
end

# ╔═╡ 12b9855f-3951-4ba9-acdf-65449b6c55bf
begin
	dgf_part(t) = scalefact_part(c,t)
	a_vals_part = dgf_part.(tgyr)
	plot(tgyr,a_vals_part,lw=2.5, label="using scalefact_part")
end

# ╔═╡ c4681d55-b75b-46c3-b5c1-789568890358
plot!(tgyr,a_vals,lw=2.5, label="using scale_fact")#They are almost identical.

# ╔═╡ dee6876a-5a84-437f-a953-881c441f34bf
plot!([13.807434539,13.807434539],[0,2.5],lc=:green,lw=2.5,label="t = 13.80743",ls=:dash)

# ╔═╡ Cell order:
# ╠═eb918f38-b7cc-11ee-1849-e9a292ae285b
# ╟─8d34d1a3-6a28-461c-9e03-cc047654ce84
# ╟─ebf05ec9-dc2e-4ac9-a231-fd05e170e0a2
# ╠═667c6bda-f78c-4923-b93c-1942f5b1f258
# ╟─dc415649-5683-415e-a1b8-a31f48080d94
# ╠═0314d0ad-d666-41fb-8e5a-bd4a0a497ba6
# ╟─a2a0d4d9-46cc-4f18-9311-de22ea571b4e
# ╠═361cb06f-4e96-48a1-a23e-e38106ef4ed8
# ╠═52930937-1f73-470d-b041-9558a11a7658
# ╟─5d11e724-326a-436f-94f4-211f99d30990
# ╠═7e4c3515-02f7-498b-829f-23e78997fb70
# ╟─3aa50f1f-04e8-4d0d-a650-c0ef8e75089f
# ╠═1473e933-d04a-4a0b-90d8-ab33df9ff3ee
# ╠═2ee6b0bf-3511-4fea-8132-9896f5f112ec
# ╟─4984ac71-b8eb-4085-a20a-9a2ac5f439a4
# ╠═2bb09fc2-f284-4d1e-8d58-55365d31dc20
# ╟─d3e187e5-13ee-4edb-bfa0-25aee07b712e
# ╠═b02a5df0-b196-4418-84b0-52633d6e49e0
# ╟─511ade9b-ccf8-4edd-a5c5-edf2ba869815
# ╠═7ed40ab7-c190-4bf0-bcf6-eccbb422784e
# ╟─ac2675fd-b0c9-4ba7-95eb-79f72bddc678
# ╠═f1b373f4-43bd-4312-bcaf-49ff3fc2d31a
# ╠═efabfaa1-f0bc-466a-bb98-49795c44e200
# ╠═b9dc315c-5a35-41fd-b27d-f1f9a4fc7e81
# ╠═1e09ba10-91d6-46f0-bc8f-b9b1a82220fd
# ╠═12b9855f-3951-4ba9-acdf-65449b6c55bf
# ╠═c4681d55-b75b-46c3-b5c1-789568890358
# ╠═dee6876a-5a84-437f-a953-881c441f34bf
