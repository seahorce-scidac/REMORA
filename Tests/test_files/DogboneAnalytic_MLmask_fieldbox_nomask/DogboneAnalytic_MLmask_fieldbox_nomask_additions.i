# The same field test without a land mask, which AMReX would otherwise evaluate on its own.
remora.max_step = 0
remora.mask_type = "none"
remora.coastbox.value_greater = -1.0
remora.coastbox.field_name    = temp
