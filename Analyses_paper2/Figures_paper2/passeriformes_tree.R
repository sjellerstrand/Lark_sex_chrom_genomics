

## Export variables and load libraries
rm(list=ls())
library(ape)
library(phytools)

tree <- read.nexus("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results/Figure-1.tre")
plot(tree)

# Drop Outgroup
tree <- drop.tip(tree, which(tree$tip.label == "Cariama_cristata"))
tree <- drop.tip(tree, which(tree$tip.label == "Micrastur_semitorquatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Caracara_cheriway"))
tree <- drop.tip(tree, which(tree$tip.label == "Microhierax_erythrogenys"))
tree <- drop.tip(tree, which(tree$tip.label == "Falco_peregrinus"))
tree <- drop.tip(tree, which(tree$tip.label == "Strigops_habroptila"))
tree <- drop.tip(tree, which(tree$tip.label == "Nestor_notabilis"))
tree <- drop.tip(tree, which(tree$tip.label == "Cacatua_sulphurea"))
tree <- drop.tip(tree, which(tree$tip.label == "Psittacula_alexandri"))
tree <- drop.tip(tree, which(tree$tip.label == "Melopsittacus_undulatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Amazona_autumnalis"))
tree <- drop.tip(tree, which(tree$tip.label == "Psittacus_erithacus"))
tree <- drop.tip(tree, which(tree$tip.label == "Xenicus_gilviventris"))

# Acanthsitti

# Tyranni
tree$tip.label[which(tree$tip.label == "Eurylaimus_ochromalus")] <- "Serilophus_lunatus"
tree$tip.label[which(tree$tip.label == "Dendrocolaptes_sanctithomae")] <- "Xiphorhynchus_elegans"
tree$tip.label[which(tree$tip.label == "Thamnophilus_doliatus")] <- "Willisornis_vidua_nigrigula"
tree$tip.label[which(tree$tip.label == "Pipra_filicauda")] <- "Chiroxiphia_lanceolata"
tree$tip.label[which(tree$tip.label == "Tyrannus_albogularis")] <- "Pyrocephalus_nanus"
tree <- drop.tip(tree, which(tree$tip.label == "Pitta_cyanea"))
tree <- drop.tip(tree, which(tree$tip.label == "Pitta_erythrogaster"))
tree <- drop.tip(tree, which(tree$tip.label == "Calyptomena_viridis"))
tree <- drop.tip(tree, which(tree$tip.label == "Smithornis_rufolateralis"))
tree <- drop.tip(tree, which(tree$tip.label == "Philepitta_castanea"))
tree <- drop.tip(tree, which(tree$tip.label == "Xenops_minutus"))
tree <- drop.tip(tree, which(tree$tip.label == "Furnarius_rufus"))
tree <- drop.tip(tree, which(tree$tip.label == "Sclerurus_rufigularis"))
tree <- drop.tip(tree, which(tree$tip.label == "Rhinocrypta_lanceolata"))
tree <- drop.tip(tree, which(tree$tip.label == "Grallaria_haplonota"))
tree <- drop.tip(tree, which(tree$tip.label == "Pygiptila_stellaris"))
tree <- drop.tip(tree, which(tree$tip.label == "Conopophaga_aurita"))
tree <- drop.tip(tree, which(tree$tip.label == "Melanopareia_maranonica"))
tree <- drop.tip(tree, which(tree$tip.label == "Cotinga_nattererii"))
tree <- drop.tip(tree, which(tree$tip.label == "Schiffornis_major"))
tree <- drop.tip(tree, which(tree$tip.label == "Tityra_cayana"))
tree <- drop.tip(tree, which(tree$tip.label == "Oxyruncus_cristatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Onychorhynchus_coronatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Myiobius_villosus"))
tree <- drop.tip(tree, which(tree$tip.label == "Piprites_chloris"))
tree <- drop.tip(tree, which(tree$tip.label == "Platyrinchus_saturatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Tachuris_rubrigastra"))
tree <- drop.tip(tree, which(tree$tip.label == "Nephelomyias_lintoni"))
tree <- drop.tip(tree, which(tree$tip.label == "Mionectes_striaticollis"))
tree <- drop.tip(tree, which(tree$tip.label == "Poecilotriccus_latirostris"))
tree <- drop.tip(tree, which(tree$tip.label == "Neopipo_cinnamomea"))

# Outer Passeri
tree$tip.label[which(tree$tip.label == "Climacteris_melanurus")] <- "Climacteris_rufus"
tree$tip.label[which(tree$tip.label == "Malurus_alboscapulatus")] <- "Malurus_cyaneus_samueli"
tree$tip.label[which(tree$tip.label == "Acanthiza_cinerea")] <- "Acanthiza_pusilla"
tree$tip.label[which(tree$tip.label == "Meliphaga_montana")] <- "Lichenostomus_cassidix"
tree$tip.label[which(tree$tip.label == "Orthonyx_temminckii")] <- "Orthonyx_spaldingii"
tree <- drop.tip(tree, which(tree$tip.label == "Atrichornis_rufescens"))
tree <- drop.tip(tree, which(tree$tip.label == "Ptilonorhynchus_violaceus"))
tree <- drop.tip(tree, which(tree$tip.label == "Ailuroedus_buccoides"))
tree <- drop.tip(tree, which(tree$tip.label == "Cormobates_leucophaea"))
tree <- drop.tip(tree, which(tree$tip.label == "Dasyornis_broadbenti"))
tree <- drop.tip(tree, which(tree$tip.label == "Pardalotus_striatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Acanthorhynchus_tenuirostris"))
tree <- drop.tip(tree, which(tree$tip.label == "Foulehaio_carunculatus"))

# Corvides
tree$tip.label[which(tree$tip.label == "Corvus_corax")] <- "Coloeus_monedula"
tree$tip.label[which(tree$tip.label == "Rhipidura_javanica")] <- "Rhipidura_dahli"
tree <- drop.tip(tree, which(tree$tip.label == "Lanius_excubitor"))
tree <- drop.tip(tree, which(tree$tip.label == "Platylophus_galericulatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Pyrrhocorax_pyrrhocorax"))
tree <- drop.tip(tree, which(tree$tip.label == "Melampitta_lugubris"))
tree <- drop.tip(tree, which(tree$tip.label == "Struthidea_cinerea"))
tree <- drop.tip(tree, which(tree$tip.label == "Phonygammus_keraudrenii"))
tree <- drop.tip(tree, which(tree$tip.label == "Paradisaea_minor"))
tree <- drop.tip(tree, which(tree$tip.label == "Metabolus_takatsukasae"))
tree <- drop.tip(tree, which(tree$tip.label == "Hypothymis_azurea"))#
tree <- drop.tip(tree, which(tree$tip.label == "Dicrurus_aeneus"))
tree <- drop.tip(tree, which(tree$tip.label == "Lamprolia_victoriae"))
tree <- drop.tip(tree, which(tree$tip.label == "Prionops_plumatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Tephrodornis_virgatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Mystacornis_crossleyi"))
tree <- drop.tip(tree, which(tree$tip.label == "Dyaphorophyia_castanea"))
tree <- drop.tip(tree, which(tree$tip.label == "Pityriasis_gymnocephala"))
tree <- drop.tip(tree, which(tree$tip.label == "Aegithina_lafresnayei"))
tree <- drop.tip(tree, which(tree$tip.label == "Tchagra_senegalus"))
tree <- drop.tip(tree, which(tree$tip.label == "Melloria_quoyi"))
tree <- drop.tip(tree, which(tree$tip.label == "Strepera_graculina"))
tree <- drop.tip(tree, which(tree$tip.label == "Peltops_blainvillii"))
tree <- drop.tip(tree, which(tree$tip.label == "Artamus_cinereus"))
tree <- drop.tip(tree, which(tree$tip.label == "Machaerirhynchus_nigripectus"))
tree <- drop.tip(tree, which(tree$tip.label == "Vireo_solitarius"))
tree <- drop.tip(tree, which(tree$tip.label == "Pachycephala_vitiensis"))
tree <- drop.tip(tree, which(tree$tip.label == "Erpornis_zantholeuca"))
tree <- drop.tip(tree, which(tree$tip.label == "Pteruthius_aeralatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Colluricincla_harmonica"))
tree <- drop.tip(tree, which(tree$tip.label == "Oriolus_chinensis"))
tree <- drop.tip(tree, which(tree$tip.label == "Pitohui_dichrous"))
tree <- drop.tip(tree, which(tree$tip.label == "Oreoica_gutturalis"))
tree <- drop.tip(tree, which(tree$tip.label == "Psophodes_cristatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Mohoua_albicilla"))
tree <- drop.tip(tree, which(tree$tip.label == "Pericrocotus_divaricatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Edolisoma_tenuirostre"))
tree <- drop.tip(tree, which(tree$tip.label == "Cinclosoma_punctatum"))
tree <- drop.tip(tree, which(tree$tip.label == "Ptilorrhoa_leucosticta"))
tree <- drop.tip(tree, which(tree$tip.label == "Pomatostomus_superciliosus"))

# Passerides Outer
tree$tip.label[which(tree$tip.label == "Petroica_multicolor")] <- "Eopsaltria_australis"
tree <- drop.tip(tree, which(tree$tip.label == "Cnemophilus_loriae"))
tree <- drop.tip(tree, which(tree$tip.label == "Loboparadisea_sericea"))
tree <- drop.tip(tree, which(tree$tip.label == "Melanocharis_versteri"))
tree <- drop.tip(tree, which(tree$tip.label == "Toxorhamphus_novaeguineae"))
tree <- drop.tip(tree, which(tree$tip.label == "Oedistoma_iliolophus"))
tree <- drop.tip(tree, which(tree$tip.label == "Notiomystis_cincta"))
tree <- drop.tip(tree, which(tree$tip.label == "Devioeca_papuana"))
tree <- drop.tip(tree, which(tree$tip.label == "Eupetes_macrocerus"))
tree <- drop.tip(tree, which(tree$tip.label == "Picathartes_gymnocephalus"))
tree <- drop.tip(tree, which(tree$tip.label == "Chaetops_frenatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Philesturnus_carunculatus"))

# Muscicapida
tree$tip.label[which(tree$tip.label == "Sitta_europea")] <- "Sitta_carolinensis"
tree$tip.label[which(tree$tip.label == "Cinclus_pallasii")] <- "Cinclus_cinclus"
tree$tip.label[which(tree$tip.label == "Sturnus_vulgaris")] <- "Lamprotornis_superbus"
tree$tip.label[which(tree$tip.label == "Muscicapa_striata")] <- "Luscinia_luscinia"
tree <- bind.tip(tree, "Ficedula_albicollis", 14.19, which(tree$tip.label=="Luscinia_luscinia"), 14.19)
tree$tip.label[which(tree$tip.label == "Turdus_albicollis")] <- "Catharus_ustulatus"
tree <- drop.tip(tree, which(tree$tip.label == "Dulus_dominicus"))
tree <- drop.tip(tree, which(tree$tip.label == "Bombycilla_garrulus"))
tree <- drop.tip(tree, which(tree$tip.label == "Ptilogonys_caudatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Hylocitrea_bonensis"))
tree <- drop.tip(tree, which(tree$tip.label == "Moho_nobilis"))
tree <- drop.tip(tree, which(tree$tip.label == "Hypocolius_ampelinus"))
tree <- drop.tip(tree, which(tree$tip.label == "Elachura_formosa"))
tree <- drop.tip(tree, which(tree$tip.label == "Mimus_polyglottos"))
tree <- drop.tip(tree, which(tree$tip.label == "Regulus_regulus"))
tree <- drop.tip(tree, which(tree$tip.label == "Regulus_calendula"))
tree <- drop.tip(tree, which(tree$tip.label == "Tichodroma_muraria"))
tree <- drop.tip(tree, which(tree$tip.label == "Troglodytes_troglodytes"))
tree <- drop.tip(tree, which(tree$tip.label == "Polioptila_caerulea"))
tree <- drop.tip(tree, which(tree$tip.label == "Salpornis_spilonota"))

# Passerida
tree$tip.label[which(tree$tip.label == "Leptocoma_sperata")] <- "Leptocoma_aspasia"
tree$tip.label[which(tree$tip.label == "Vidua_regia")] <- "Vidua_chalybeata"
tree$tip.label[which(tree$tip.label == "Prunella_fulvescens")] <- "Prunella_stophiata"
tree$tip.label[which(tree$tip.label == "Motacilla_alba")] <- "Motacilla_alba_alba"
tree$tip.label[which(tree$tip.label == "Fringilla_montifringilla")] <- "Haemorhous_mexicanus"
tree$tip.label[which(tree$tip.label == "Emberzia_citrinella")] <- "Emberzia_elegans"
tree$tip.label[which(tree$tip.label == "Coereba_flaveola")] <- "Camarhynchus_parvulus"
tree$tip.label[which(tree$tip.label == "Passerella_iliaca")] <- "Passeculus_sandwichensis"
tree$tip.label[which(tree$tip.label == "Sturnella_neglecta")] <- "Molothrus_aeneus"
tree$tip.label[which(tree$tip.label == "Seiurus_aurocapilla")] <- "Geothlypis trichas"
tree$tip.label[which(tree$tip.label == "Estrilda_melpoda")] <- "Taeniopygia_guttata"
tree <- drop.tip(tree, which(tree$tip.label == "Promerops_gurneyi"))
tree <- drop.tip(tree, which(tree$tip.label == "Modulatrix_stictigula"))
tree <- drop.tip(tree, which(tree$tip.label == "Dicaeum_hypoleucum"))
tree <- drop.tip(tree, which(tree$tip.label == "Irena_puella"))
tree <- drop.tip(tree, which(tree$tip.label == "Irena_cyanogastra"))
tree <- drop.tip(tree, which(tree$tip.label == "Chloropsis_cochinchinensis"))
tree <- drop.tip(tree, which(tree$tip.label == "Chloropsis_sonnerati"))
tree <- drop.tip(tree, which(tree$tip.label == "Peucedramus_taeniatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Urocynchramus_pylzowi"))
tree <- drop.tip(tree, which(tree$tip.label == "Ploceus_cucullatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Calcarius_lapponicus"))
tree <- drop.tip(tree, which(tree$tip.label == "Mitrospingus_cassinii"))
tree <- drop.tip(tree, which(tree$tip.label == "Tangara_arthus"))
tree <- drop.tip(tree, which(tree$tip.label == "Spizella_passerina"))
tree <- drop.tip(tree, which(tree$tip.label == "Ammodramus_savannarum"))
tree <- drop.tip(tree, which(tree$tip.label == "Icteria_virens"))
tree <- drop.tip(tree, which(tree$tip.label == "Icterus_cucullatus"))
tree <- drop.tip(tree, which(tree$tip.label == "Setophaga_magnolia"))
tree <- drop.tip(tree, which(tree$tip.label == "Calyptophilus_tertius"))
tree <- drop.tip(tree, which(tree$tip.label == "Zeledonia_coronata"))
tree <- drop.tip(tree, which(tree$tip.label == "Phaenicophilus_palmarum"))
tree <- drop.tip(tree, which(tree$tip.label == "Microligea_palustris"))
tree <- drop.tip(tree, which(tree$tip.label == "Spindalis_zena"))
tree <- drop.tip(tree, which(tree$tip.label == "Nesospingus_speculiferus"))

# Passerides Outer
tree$tip.label[which(tree$tip.label == "Remiz_consobrinus")] <- "Remiz_coronatus"
tree <- bind.tip(tree, "Sylviparus_modestus", 15.22124, which(tree$tip.label=="Parus_major"), 15.22124)
tree <- bind.tip(tree, "Melanochlora_sultanea", 13.85183, which(tree$tip.label=="Parus_major"), 13.85183)
tree <- bind.tip(tree, "Poecile_montanus", 10.77034, which(tree$tip.label=="Parus_major"), 10.77034)
tree <- drop.tip(tree, which(tree$tip.label == "Hyliota_flavigaster"))
tree <- drop.tip(tree, which(tree$tip.label == "Culicicapa_ceylonensis"))
tree <- drop.tip(tree, which(tree$tip.label == "Chelidorhynx_hypoxanthus"))

# Sylvioidea
tree$tip.label[which(tree$tip.label == "Eremophila_alpestris")] <- "Alauda_arvensis"
tree <- drop.tip(tree, which(tree$tip.label == "Panurus_biarmicus"))
tree <- drop.tip(tree, which(tree$tip.label == "Nicator_chloris"))
tree <- drop.tip(tree, which(tree$tip.label == "Macrosphenus_flavicans"))
tree <- drop.tip(tree, which(tree$tip.label == "Sylvietta_virens"))
tree <- drop.tip(tree, which(tree$tip.label == "Neomixis_viridis"))
tree <- drop.tip(tree, which(tree$tip.label == "Orthotomus_castaneiceps"))
tree <- drop.tip(tree, which(tree$tip.label == "Cisticola_anonymus"))
tree <- drop.tip(tree, which(tree$tip.label == "Graueria_vittata"))
tree <- drop.tip(tree, which(tree$tip.label == "Nesillas_typica"))
tree <- drop.tip(tree, which(tree$tip.label == "Acrocephalus_orientalis"))
tree <- drop.tip(tree, which(tree$tip.label == "Donacobius_atricapilla"))
tree <- drop.tip(tree, which(tree$tip.label == "Oxylabes_madagascariensis"))
tree <- drop.tip(tree, which(tree$tip.label == "Xanthomixis_zosterops"))
tree <- drop.tip(tree, which(tree$tip.label == "Robsonius_rabori"))
tree <- drop.tip(tree, which(tree$tip.label == "Locustella_lanceolata"))
tree <- drop.tip(tree, which(tree$tip.label == "Megalurus_palustris"))
tree <- drop.tip(tree, which(tree$tip.label == "Pnoepyga_pusilla"))
tree <- drop.tip(tree, which(tree$tip.label == "Hirundo_rustica"))
tree <- drop.tip(tree, which(tree$tip.label == "Progne_subis"))
tree <- drop.tip(tree, which(tree$tip.label == "Bleda_syndactylus"))
tree <- drop.tip(tree, which(tree$tip.label == "Curruca_nana"))
tree <- drop.tip(tree, which(tree$tip.label == "Cholornis_unicolor"))
tree <- drop.tip(tree, which(tree$tip.label == "Zosterops_everetti"))
tree <- drop.tip(tree, which(tree$tip.label == "Timalia_pileata"))
tree <- drop.tip(tree, which(tree$tip.label == "Trichastoma_pyrrogenys"))
tree <- drop.tip(tree, which(tree$tip.label == "Leiothrix_lutea"))
tree <- drop.tip(tree, which(tree$tip.label == "Phylloscopus_trochilus"))
tree <- drop.tip(tree, which(tree$tip.label == "Seicercus_montis"))
tree <- drop.tip(tree, which(tree$tip.label == "Pholidornis_rushiae"))
tree <- drop.tip(tree, which(tree$tip.label == "Hylia_prasina"))
tree <- drop.tip(tree, which(tree$tip.label == "Aegithalos_concinnus"))
tree <- drop.tip(tree, which(tree$tip.label == "Psaltriparus_minimus"))
tree <- drop.tip(tree, which(tree$tip.label == "Erythrocercus_mccallii"))
tree <- drop.tip(tree, which(tree$tip.label == "Pycnonotus_jocosus"))
tree <- drop.tip(tree, which(tree$tip.label == "Scotocerca_inquieta"))
tree <- drop.tip(tree, which(tree$tip.label == "Tesia_cyaniventer"))
tree <- drop.tip(tree, which(tree$tip.label == "Abroscopus_albogularis"))

plot(tree)

png("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results/Figures/Supplementary/GERP_tree.png", width=4000, height=3000, res=300)
plot(tree)
dev.off()

write.tree(tree, file = "C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results/Passeriformes_tree_GERP.nwk", append = FALSE,
           digits = 100, tree.names = FALSE)
